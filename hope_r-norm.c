#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <semaphore.h>
#include <apop.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_multifit.h>

/********************************HOPE: Sum of norm of residuals fitting for Hough Mode Extension Data*****************************************/
//Copyright (C) 2012 Daniel Strano

/*'Hope' is an open-source project based on Oberheide and Forbes GEOPHYSICAL RESEARCH LETTERS, VOL. 35, L04816, doi:10.1029/2007GL032397, 2008. The program was originally designed for the purpose of handling the challenge of empirically fitting multi-component complex wave data, in situations where the goodness of fit of a least squares model is questionable. The program was originally designed to work with three self-consistent variables to fit simultaneously. The variances of the variables are first scaled to satisify the assumptions of constant variance and independence of errors, necessary for a consistent least squares model. Then, the program takes advantage of the finite phase domain to indirectly search for a maximization of the sum of absolute (norms of) residuals.*/

/*********************************************************************************************************************************************/

#define depVar 3
const int rawblock = 29;
const int block = 33;
const int heightsOut = 8;
const int heightsIn = 8;
const int HME_N = 2;
const int HMEBlockCount = 8;
const int satrows = 264;

const double PI = 3.141592653589793;
const double PHASE_INTERVAL = 0.001;
const int PHASE_INT_INV = 1000;
const double P_I_x_PI = 0.001 * 3.141592653589793;

int outrows=0;
double totalSumSqr;
double TotalSN[depVar];

int global_phase_iter1 = 0;
int global_phase_iter2 = 0;

int threadtracker[4];

pthread_mutex_t iter_update;
pthread_mutex_t compare;

sem_t threadsAvail;

gsl_matrix *Full_Descrip;
gsl_vector *dependent;
gsl_vector *BestModel;

int minIndepOffset;
double *threadIndep;

double BestPhase1;
double BestPhase2;

double BestResid = 0;
double BestRNorm = 0;

int itemscount;

void *threadParallel(void *arg) {

	int *threadNum = (int *)arg; //Argument is used for thread scheduling (used at bottom of thread function)

	int phase_iter1, phase_iter2; //Local phase iteration holder variables

	//Lock the current global phase iteration counters, pull the next value for the local iteration, and update global counters:
	pthread_mutex_lock(&iter_update);
		phase_iter1 = global_phase_iter1;
		phase_iter2 = global_phase_iter2;
		++global_phase_iter1;
		if (global_phase_iter1 >= PHASE_INT_INV) {
			global_phase_iter1 = 0;
			++global_phase_iter2;
		}
	pthread_mutex_unlock(&iter_update);
	
	int loop1, loop2, loop3, loop4, loop5;
	int returncode;
	float phase_shift[2]= {P_I_x_PI * phase_iter1, P_I_x_PI * phase_iter2}; //Phase shifts calculated from phase iteration counters
	float Rotation[8]; //Used to rotate complex vectors

	int holder, holder2, holder3, holder4, holder5, holder6; //For loop optimization

	double chiSqr; //Sum of squares of residuals of each iterative model
	float coeff_of_determ; //Normalized coefficient of determination ""

	//Calculate rotation factors to apply phase shift corresponding to the current iteration:	
	for (loop1 = 0; loop1 < HME_N; loop1++) {
		holder = loop1 * 4;
		holder2 = holder + 2;
		for (loop2 = 0; loop2 < 2; loop2++) {
			Rotation[holder + loop2] = cos(phase_shift[loop1]);
			if (loop2) {
				Rotation[holder2 + loop2] = -1 * sin(phase_shift[loop1]);
			}
			else {
				Rotation[holder2 + loop2] = sin(phase_shift[loop1]);
			}
		}
	}

	//Declare objects for GSL least squares model for fixed phases:
	gsl_matrix *Local_Descrip = gsl_matrix_alloc(itemscount, HME_N); //Local description of amplitude-only model
	gsl_matrix *covariance = gsl_matrix_alloc(HME_N, HME_N);
	gsl_multifit_linear_workspace *modelspace = gsl_multifit_linear_alloc(itemscount, HME_N);
	gsl_vector *beta = gsl_vector_alloc(HME_N); //fitting parameters

	//Apply rotations and initialize 'Local_Descrip'
	for (loop5 = 0; loop5 < HME_N; loop5++) {
		holder = loop5 * 4;
		holder2 = holder + 1;
		holder3 = holder2 + 1;
		holder4 = holder3 + 1;
		holder6 = loop5 * minIndepOffset;
		for (loop1 = 0; loop1 < itemscount; loop1 += 2) {
			holder5 = loop1 + 1;
			gsl_matrix_set(Local_Descrip, loop1, loop5,
				Rotation[holder] * threadIndep[holder6 + loop1] + Rotation[holder4] * threadIndep[holder6 + holder5]);
			gsl_matrix_set(Local_Descrip, holder5, loop5,
				Rotation[holder2] * threadIndep[holder6 + holder5] + Rotation[holder3] * threadIndep[holder6 + loop1]);
		}
	}

	//Find least squares fit:
	returncode = gsl_multifit_linear(Local_Descrip, dependent, beta, covariance, &chiSqr, modelspace);
	coeff_of_determ = 1 - chiSqr/totalSumSqr;

	//Free model space for fit (covariance is also not used):
	gsl_matrix_free(covariance);
	gsl_multifit_linear_free(modelspace);

	//From partial model beta parameters, calculate equivalent full model parameters:
	double sumNormResid[3] = {0, 0, 0};
	double R_norm = 0;
	double R_resid, I_resid;
	double Full_Beta[2 * HME_N];
	for (loop1 = 0; loop1 < HME_N; loop1++) {
		Full_Beta[2 * loop1] = gsl_vector_get(beta, loop1) * cos(phase_shift[loop1]);
		Full_Beta[2 * loop1 + 1] = gsl_vector_get(beta, loop1) * sin(phase_shift[loop1]);
	}

	//Re-weight variables to equal total weight for sum of norms of residuals model, 
	//calculate sum of norms of residuals for full model from least squares beta coefficients:
	int maxiter1 = 2 * outrows;
	int maxiter2 = maxiter1 * depVar;
	int maxiter3 = HME_N * 2;
	for (loop2 = 0; loop2 < maxiter2; loop2 += maxiter1) {
		holder3 = (int)(loop2/maxiter1);
		for (loop1 = 0; loop1 < maxiter1; loop1 += 2) {
			R_resid = gsl_vector_get(dependent, loop1 + loop2);
			holder = loop2 + loop1;
			holder2 = holder + 1;
			for (loop3 = 0; loop3 < maxiter3; loop3++) R_resid -= Full_Beta[loop3] * gsl_matrix_get(Full_Descrip, holder, loop3);
			I_resid = gsl_vector_get(dependent, holder2);
			for (loop3 = 0; loop3 < maxiter3; loop3++) I_resid -= Full_Beta[loop3] * gsl_matrix_get(Full_Descrip, holder2, loop3);
			sumNormResid[holder3] += sqrt((R_resid * R_resid) + (I_resid * I_resid));
		}
	}

	//Calculated normalized coefficient of determination for sum of norms of residuals:
	for (loop1 = 0; loop1 < depVar; loop1++) R_norm -= sumNormResid[loop1]/TotalSN[loop1];
	R_norm = R_norm/3;
	R_norm++;


	//Lock global best model comparisons and updates,
	//compare this thread's "R-norm" and update global best model if a better model has been found:
	pthread_mutex_lock(&compare);
		if (BestRNorm < R_norm) {
			BestRNorm = R_norm;
			BestResid = coeff_of_determ;
			BestPhase1 = phase_shift[0];
			BestPhase2 = phase_shift[1];
			for (loop1 = 0; loop1 < HME_N; loop1++) {			
				gsl_vector_set(BestModel, loop1, gsl_vector_get(beta, loop1));
			}
		}
	pthread_mutex_unlock(&compare);

	//Free remaining manually allocated memory objects:
	gsl_matrix_free(Local_Descrip);
	gsl_vector_free(beta);

	//Prepare to signal to loop in main() that this thread identity is becoming free, then post to semaphore:
	threadtracker[*threadNum] = 0;
	sem_post(&threadsAvail);

	return NULL;
}

int main(int argc, char *argv[])
{

	//Thread expects input files as specified below:
	if (argc != 4) {
		printf("usage: hope [saber_t] [tidi_e] [tidi_n]\n");
		printf("Note: Satellite data should be in *.csv form, in a single column, aligned by row, with no column or row names. Dummy value is -999. \n");
	}
	else {

		gsl_matrix *satellite[3];

		double *minIndep[2];

		int lev_rows;
		double height[heightsIn];
		double latitudes[rawblock];

		double *HME_ptr[8][HME_N];

		double **uninterp[HME_N][depVar*2];
		double **interp_lat[HME_N][depVar*2];
		double **interp_lev[HME_N][depVar*2];

		BestModel = gsl_vector_alloc(2);

		pthread_t zero;
		pthread_t one;
		pthread_t two;
		pthread_t three;

		sem_init(&threadsAvail, 0, 4);

		long long maxpage;
		long long pagein = 1;
		int comp, comp_int, tenth_percent_complete = 0, percent_complete = 0;
	
		int loop, loop1, loop2, loop3, loop4, loop5, loop6;

		int threadID[4];

		double TotalSumSqr[depVar];
		double Variance[depVar];
		double AbsDeviance[depVar];
		double Scale[depVar], DevScale[depVar];

		double scaledRow;

		int holder;
		int holder2;
		 
		for (loop = 0; loop < depVar; loop++) {
			TotalSumSqr[loop] = 0;
			TotalSN[loop] = 0;
			Variance[loop] = 0;  //optional
			Scale[loop] = 0; //optional
		}

		for (loop = 0; loop < 4; loop++) {
			threadtracker[loop] = 0;
			threadID[loop] = loop;
		}

			//Read in raw data:
	
				//(Independent) descriptor data:
		apop_text_to_db(.text_file="HME1.csv", .tabname="X1", .has_row_names = 0, .has_col_names = 1);
		HME_ptr[0][0] = apop_vector_to_array(apop_query_to_vector("SELECT lat FROM X1"));
		HME_ptr[1][0] = apop_vector_to_array(apop_query_to_vector("SELECT lev FROM X1"));
		HME_ptr[2][0] = apop_vector_to_array(apop_query_to_vector("SELECT amp_t FROM X1"));
		HME_ptr[3][0] = apop_vector_to_array(apop_query_to_vector("SELECT phase_t FROM X1"));
		HME_ptr[4][0] = apop_vector_to_array(apop_query_to_vector("SELECT amp_e FROM X1"));
		HME_ptr[5][0] = apop_vector_to_array(apop_query_to_vector("SELECT phase_e FROM X1"));
		HME_ptr[6][0] = apop_vector_to_array(apop_query_to_vector("SELECT amp_n FROM X1"));
		HME_ptr[7][0] = apop_vector_to_array(apop_query_to_vector("SELECT phase_n FROM X1"));
	
		apop_text_to_db(.text_file="HME2.csv", .tabname="X2", .has_row_names = 0, .has_col_names = 1);
		HME_ptr[0][1] = apop_vector_to_array(apop_query_to_vector("SELECT lat FROM X2"));
		HME_ptr[1][1] = apop_vector_to_array(apop_query_to_vector("SELECT lev FROM X2"));
		HME_ptr[2][1] = apop_vector_to_array(apop_query_to_vector("SELECT amp_t FROM X2"));
		HME_ptr[3][1] = apop_vector_to_array(apop_query_to_vector("SELECT phase_t FROM X2"));
		HME_ptr[4][1] = apop_vector_to_array(apop_query_to_vector("SELECT amp_e FROM X2"));
		HME_ptr[5][1] = apop_vector_to_array(apop_query_to_vector("SELECT phase_e FROM X2"));
		HME_ptr[6][1] = apop_vector_to_array(apop_query_to_vector("SELECT amp_n FROM X2"));
		HME_ptr[7][1] = apop_vector_to_array(apop_query_to_vector("SELECT phase_n FROM X2"));

				//(Dependent) variable data:	
		apop_text_to_db(.text_file=argv[1], .tabname="Y_t", .has_row_names = 0, .has_col_names = 1);
		satellite[0] = apop_query_to_matrix("SELECT * FROM Y_t");
	
		apop_text_to_db(.text_file=argv[2], .tabname="Y_e", .has_row_names = 0, .has_col_names = 1);
		satellite[1] = apop_query_to_matrix("SELECT * FROM Y_e");

		apop_text_to_db(.text_file=argv[3], .tabname="Y_n", .has_row_names = 0, .has_col_names = 1);
		satellite[2] = apop_query_to_matrix("SELECT * FROM Y_n");


		//Need to know how many non-dummy rows are in satellite data
		//Calculate total sum of squares and total sum of norms of satellite data while tallying non-dummy rows in 'outrows':
		outrows = 0;
		for (loop1 = 0; loop1 < satrows; loop1++) {
			if (gsl_matrix_get(satellite[0], loop1, 0) > -990) {
				if (gsl_matrix_get(satellite[1], loop1, 0) > -990) {
					if (gsl_matrix_get(satellite[2], loop1, 0) > -990) {
						for (loop2 = 0; loop2 < depVar; loop2++) {
							for (loop3 = 0; loop3 < 2; loop3++) {
								TotalSumSqr[loop2] += (gsl_matrix_get(satellite[loop2], loop1, loop3) * gsl_matrix_get(satellite[loop2], loop1, loop3)); 
							}
							TotalSN[loop2] += sqrt((gsl_matrix_get(satellite[loop2], loop1, 0) * gsl_matrix_get(satellite[loop2], loop1, 0)) + (gsl_matrix_get(satellite[loop2], loop1, 1) * gsl_matrix_get(satellite[loop2], loop1, 1)));
						}
						++outrows;
					}
				}
			}
		}

		itemscount = outrows * depVar * 2; //Total number of data items (rows * columns)

		//Allocate arrays for uninterpolated data, lattitude interpolated data, and data interpolated by both lattitude and elevation:
		for (loop1 = 0; loop1 < HME_N; loop1++) {

			for (loop2 = 0; loop2 < depVar * 2; loop2++) {
				uninterp[loop1][loop2]   = malloc(HMEBlockCount * sizeof(double *));
				interp_lat[loop1][loop2] = malloc(block * sizeof(double *));
				interp_lev[loop1][loop2] = malloc(block * sizeof(double *));
				for (loop3 = 0; loop3 < HMEBlockCount; loop3++) {
					uninterp[loop1][loop2][loop3]   = malloc(rawblock * sizeof(double));
				}
				for (loop3 = 0; loop3 < block; loop3++) {
					interp_lat[loop1][loop2][loop3] = malloc(HMEBlockCount * sizeof(double));
					interp_lev[loop1][loop2][loop3] = malloc(heightsOut * sizeof(double));
				}
			}
		}

		//Initialize block of uninterpolated data:
		for (loop1 = 0; loop1 < HMEBlockCount; loop1++) {
			holder = rawblock * loop1;
			for (loop2 = 0; loop2 < 2 * depVar; loop2 += 2) {
				for (loop3 = 0; loop3 < HME_N; loop3++) {
					for (loop4 = 0; loop4 < rawblock; loop4++) {
						uninterp[loop3][0 + loop2][loop1][(rawblock - 1) - loop4] = *(HME_ptr[2 + loop2][loop3] + (holder + loop4)) * cos(*(HME_ptr[3 + loop2][loop3] + holder + loop4) / 12 * PI);
						uninterp[loop3][1 + loop2][loop1][(rawblock - 1) - loop4] = *(HME_ptr[2 + loop2][loop3] + (holder + loop4)) * sin(*(HME_ptr[3 + loop2][loop3] + holder + loop4) / 12 * PI);
					}
				}
			}
		}

		//Read lattitudes from input:
		for (loop1 = 0; loop1 < rawblock; loop1++) {
			latitudes[(rawblock - 1) - loop1] = *(HME_ptr[0][0] + loop1);
		}
		
		//Read elevations from input:
		for (loop1 = 0; loop1 < heightsIn; loop1++) {
			height[loop1] = *(HME_ptr[1][0] + (loop1 * rawblock));
		}
		

	//GSL Spline Interpolation:

		//Lattitude interpolation:
		for (loop = 0; loop < HME_N; loop++) {		
			for (loop3 = 0; loop3 < 2; loop3++) {
				for (loop2 = 0; loop2 < HMEBlockCount; loop2++) {
					for (loop4 = 0; loop4 < depVar * 2; loop4 += 2) {
						gsl_interp_accel *lat_acc = gsl_interp_accel_alloc();
						gsl_spline *lat_spline = gsl_spline_alloc(gsl_interp_cspline, rawblock);
						gsl_spline_init(lat_spline, latitudes, uninterp[loop][loop4 + loop3][loop2], rawblock);

						for (loop1 = 0; loop1 < block; loop1++) {
							interp_lat[loop][loop4 + loop3][loop1][loop2] = gsl_spline_eval(lat_spline, (loop1 * 5) - 80, lat_acc);
						}

						gsl_spline_free (lat_spline);
						gsl_interp_accel_free (lat_acc);
					}
				}
			}
		}

		//Altitude ("lev") interpolation:
		holder = 0;
		for (loop = 0; loop < HME_N; loop++) {		
			for (loop3 = 0; loop3 < 2; loop3++) {
				for (loop1 = 0; loop1 < block; loop1++) {
					for (loop4 = 0; loop4 < depVar * 2; loop4 += 2) {
						gsl_interp_accel *lev_acc = gsl_interp_accel_alloc();
						gsl_spline *lev_spline = gsl_spline_alloc(gsl_interp_cspline, heightsIn);
						gsl_spline_init(lev_spline, height, interp_lat[loop][loop4 + loop3][loop1], heightsIn);

						for (loop2 = 0; loop2 < heightsOut; loop2++) {
							interp_lev[loop][loop4 + loop3][loop1][loop2] = gsl_spline_eval(lev_spline, 87.5 + (loop2 * 2.5), lev_acc);

						}

						gsl_spline_free (lev_spline);
						gsl_interp_accel_free (lev_acc);
					}
				}
			}
		}

		//Variance and deviance calculation (for normalizing different variables):
		for (loop1 = 0; loop1 < depVar; loop1++) {
			Variance[loop1] = TotalSumSqr[loop1] / (outrows/2);
			AbsDeviance[loop1] = TotalSN[loop1] / (outrows/2);
		}

		for (loop1 = 0; loop1 < depVar; loop1++) {
			Scale[loop1] =  sqrt(Variance[0]/ Variance[loop1]);
			DevScale[loop1] = AbsDeviance[0] / AbsDeviance[loop1];
		}

		//In interative fit, phases are held fixed while amplitudes are fit
		//Phases are incremented after each amplitude fit
		//Create an array with the minimum information from the independent matrix necessary to perform iterative fit:
		for (loop1 = 0; loop1 < HME_N; loop1++) {
			minIndep[loop1] = malloc(satrows * depVar * 2 * sizeof(double));
		}

		loop4 = 0;
		for (loop1 = 0; loop1 < heightsOut; loop1++) {
			for (loop3 = 0; loop3 < depVar; loop3++) {
				for (loop4 = 0; loop4 < block; loop4++) {
					for (loop5 = 0; loop5 < HME_N; loop5++) {
						for (loop = 0; loop < 2; loop++) {
							minIndep[loop5][(loop3 * 2 * satrows) + (block * 2 * loop1) + (loop4 * 2) + loop] = interp_lev[loop5][2 * loop3 + loop][block - 1 - loop4][loop1];
						}
					}
				}
			}
		}

		//Garbage collection after pre-processing:
		for (loop1 = 0; loop1 < HME_N; loop1++) {
			for (loop2 = 0; loop2 < depVar * 2; loop2++) {
				for (loop3 = 0; loop3 < HMEBlockCount; loop3++) {
					free(uninterp[loop1][loop2][loop3]);
				}
				for (loop3 = 0; loop3 < block; loop3++) {
					free(interp_lat[loop1][loop2][loop3]);
					free(interp_lev[loop1][loop2][loop3]);
				}
				free(uninterp[loop1][loop2]);
				free(interp_lat[loop1][loop2]);
				free(interp_lev[loop1][loop2]);				
			}
		}

		//Dummy values in satellite data can not be passed to GSL fitting functions
		//Real and imaginary components alternate.
		//The dependent variables counted by "depVar" are concatenated as contigious blocks
		//Remove dummy values and concatenate satellite vectors;
		dependent = gsl_vector_alloc(2 * depVar * outrows);
		totalSumSqr = 0;
		for (loop2 = 0; loop2 < depVar; loop2++) {
			loop5 = 0;
			for (loop1 = 0; loop1 < satrows; loop1++) {
				if (gsl_matrix_get(satellite[0], loop1, 0) > -990) {
					if (gsl_matrix_get(satellite[1], loop1, 0) > -990) {
						if (gsl_matrix_get(satellite[2], loop1, 0) > -990) {
							for (loop3 = 0; loop3 < 2; loop3++) {
								scaledRow = Scale[loop2] * gsl_matrix_get(satellite[loop2], loop1, loop3);
								gsl_vector_set(dependent, (2 * outrows * loop2) + loop5 + loop3, scaledRow);
								totalSumSqr += scaledRow * scaledRow;
								for (loop4 = 0; loop4 < HME_N; loop4++){
									minIndep[loop4][(2 * outrows * loop2) + loop5 + loop3] = Scale[loop2] * minIndep[loop4][(2 * satrows * loop2) + 2 * loop1 + loop3];
								}
							}
							loop5+=2;
						}
					}
				}
			}
		}

		//Reallocate shrink minIndep once dummy values are removed:
		int *error_ptr;
		minIndepOffset = outrows * depVar * 2;
		for (loop1 = 0; loop1 < HME_N; loop1++) {
			 error_ptr = realloc(minIndep[loop1], sizeof(double) * minIndepOffset);
		}

		//Create a 1D array of the minimum data that needs to be passed to threads for iterative fit:
		threadIndep = malloc(sizeof(double) * minIndepOffset * HME_N);
		for (loop1 = 0; loop1 < HME_N; loop1++) {
			for (loop2 = 0; loop2 < minIndepOffset; loop2++) {
				threadIndep[minIndepOffset * loop1 + loop2] = minIndep[loop1][loop2];
			}
		}

		//Create an array with the full 'clean' independent data for an (analytical) R^2 fit to compare against:
		Full_Descrip = gsl_matrix_alloc((outrows * depVar) * 2, HME_N * 2);
		gsl_matrix *covariance = gsl_matrix_alloc(HME_N * 2, HME_N * 2);
		gsl_multifit_linear_workspace *modelspace = gsl_multifit_linear_alloc((outrows * depVar) * 2, HME_N * 2);
		gsl_vector *beta = gsl_vector_alloc(HME_N * 2);

		for (loop1 = 0; loop1 < 2 * outrows * depVar; loop1++) {
			for (loop2 = 0; loop2 < HME_N; loop2++) {
				for (loop3 = 0; loop3 < 2; loop3++) {
					if (loop3 && !(loop1 % 2)) {
						gsl_matrix_set(Full_Descrip, loop1, 2 * loop2 + loop3, -minIndep[loop2][loop1 + loop3]);
					}
					else {
						gsl_matrix_set(Full_Descrip, loop1, 2 * loop2 + loop3, minIndep[loop2][loop1 + loop3]);
					}
				}
			}
		}

		//'minIndep' is no longer needed
		//Release 'minIndep':
		for (loop1 = 0; loop1 < 2; loop1++) {
			free(minIndep[loop1]);
		}

		//Find R^2 fit for full independent array:
		double chiSqr;
		gsl_multifit_linear(Full_Descrip, dependent, beta, covariance, &chiSqr, modelspace);
		double coeff_of_determ = 1 - chiSqr/totalSumSqr;

		//Find sum of absolute norms of residuals for the returned least squares model (normalizing variables weights to new model moment):
		double sumNormResid[3] = {0, 0, 0};
		double R_norm = 0;
		double R_resid, I_resid;
		double Full_Beta[2 * HME_N];
		double amp_phase[2 * HME_N];
		for (loop1 = 0; loop1 < HME_N * 2; loop1++) {
			Full_Beta[loop1] = gsl_vector_get(beta, loop1);
		}

		for (loop1 = 0; loop1 < HME_N * 2; loop1+=2) {
			amp_phase[loop1] = sqrt(Full_Beta[loop1] * Full_Beta[loop1] + Full_Beta[loop1 + 1] * Full_Beta[loop1 + 1]);
			amp_phase[loop1+1] = atan(Full_Beta[loop1 + 1] / Full_Beta[loop1]);
		}

		int maxiter1 = 2 * outrows;
		int maxiter2 = maxiter1 * depVar;
		for (loop2 = 0; loop2 < maxiter2; loop2 += maxiter1) {
			for (loop1 = 0; loop1 < maxiter1; loop1 += 2) {
				R_resid = gsl_vector_get(dependent, loop1 + loop2);
				for (loop3 = 0; loop3 < HME_N * 2; loop3++) R_resid -= Full_Beta[loop3] * gsl_matrix_get(Full_Descrip, loop2 + loop1, loop3);
				I_resid = gsl_vector_get(dependent, loop1 + loop2 + 1);
				for (loop3 = 0; loop3 < HME_N * 2; loop3++) I_resid -= Full_Beta[loop3] * gsl_matrix_get(Full_Descrip, loop2 + loop1 + 1, loop3);
				sumNormResid[(int)(loop2/maxiter1)] += sqrt((R_resid * R_resid) + (I_resid * I_resid));
			}
		}

		for (loop1 = 0; loop1 < depVar; loop1++) R_norm -= sumNormResid[loop1]/TotalSN[loop1];
		R_norm = R_norm/3;
		R_norm += 1;

		//Pre-processing has completed
		
		//Start calling iterative fitting threads for R-norm:
		sem_wait(&threadsAvail);
		threadtracker[0] = 1;
		pthread_create(&zero,NULL,threadParallel, &threadID[0]);
		++pagein;

		sem_wait(&threadsAvail);
		threadtracker[1] = 1;
		pthread_create(&one,NULL,threadParallel, &threadID[1]);
		++pagein;

		sem_wait(&threadsAvail);
		threadtracker[2] = 1;
		pthread_create(&two,NULL,threadParallel, &threadID[2]);
		++pagein;

		sem_wait(&threadsAvail);
		threadtracker[3] = 1;
		pthread_create(&three,NULL,threadParallel, &threadID[3]);
		++pagein;

		//With every thread initially called, begin thread scheduling:

		maxpage = (PHASE_INT_INV) * (PHASE_INT_INV);
	
		comp = floor(maxpage/1000);
		comp_int = comp;
	
		while (pagein < maxpage) {

			sem_wait(&threadsAvail);
			if (!threadtracker[0]) {
				threadtracker[0] = 1;
				pthread_join(zero, NULL);
				pthread_create(&zero,NULL,threadParallel, &threadID[0]);
			}
			else if (!threadtracker[1]) {
				threadtracker[1] = 1;
				pthread_join(one, NULL);
				pthread_create(&one,NULL,threadParallel, &threadID[1]);
			}
			else if (!threadtracker[2]) {
				threadtracker[2] = 1;
				pthread_join(two, NULL);
				pthread_create(&two,NULL,threadParallel, &threadID[2]);
			}
			else {
				threadtracker[3] = 1;
				pthread_join(three, NULL);
				pthread_create(&three, NULL, threadParallel, &threadID[3]);
			}
			++pagein;

			//Output progress on percent or 1/10 of a percent of completion:
			if (pagein > comp) {
				comp = comp + comp_int;
				++tenth_percent_complete;
				if (tenth_percent_complete > 9) {
					++percent_complete;
					tenth_percent_complete = 0;	
					//printf("Phase intercept 1: %f \n", BestPhase1);
					//printf("Phase intercept 2: %f \n", BestPhase2);
					//printf("R Squared average: %f \n", BestResid);
					//printf("%d.%d percent complete \n \n", percent_complete, tenth_percent_complete);
				}
				
				printf("Amplitude 1: %f \n", gsl_vector_get(BestModel, 0));
				printf("Phase intercept 1: %f \n", BestPhase1);
				printf("Amplitude 2: %f \n", gsl_vector_get(BestModel, 1));
				printf("Phase intercept 2: %f \n", BestPhase2);
				printf("Best R norm: %f \n", BestRNorm);
				printf("Corresponding R^2: %f \n", BestResid);
				printf("Best R^2: %f \n", coeff_of_determ);
				printf("Best R^2 Amplitude 1: %f \n", amp_phase[0]);
				printf("Best R^2 Phase intercept 1: %f \n", amp_phase[1]);
				printf("Best R^2 Amplitude 2: %f \n", amp_phase[2]);
				printf("Best R^2 Phase intercept 2: %f \n", amp_phase[3]);
				printf("R norm corresponding to analytical best R^2: %f \n", R_norm);
				printf("%d.%d percent complete \n \n", percent_complete, tenth_percent_complete);
				
			}
		}

		//All iterations have been called
		//Join last calls:
		pthread_join(zero, NULL);
		pthread_join(one, NULL);
		pthread_join(two, NULL);
		pthread_join(three, NULL);

		//Print final result:
		printf("Best phase 1 intercept: %f \n", BestPhase1);
		printf("Best phase 2 intercept: %f \n", BestPhase2);
		printf("Best R norm: %f \n", BestRNorm);
		printf("Best corresponding R squared: %f \n", BestResid);
		printf("Corresponding R^2: %f \n", BestResid);
		printf("Best R^2: %f \n", coeff_of_determ);
		printf("Best R^2 Amplitude 1: %f \n", amp_phase[0]);
		printf("Best R^2 Phase intercept 1: %f \n", amp_phase[1]);
		printf("Best R^2 Amplitude 2: %f \n", amp_phase[2]);
		printf("Best R^2 Phase intercept 2: %f \n", amp_phase[3]);
		printf("R norm corresponding to analytical best R^2: %f \n", R_norm);

		//Free remaining manually allocated memory objects:
		for (loop1 = 0; loop1 < depVar; loop1++) {
			gsl_matrix_free(satellite[loop1]);
		}
		gsl_vector_free(dependent);
		gsl_vector_free(BestModel);
		gsl_matrix_free(Full_Descrip);
		
	}
	return 0;
}
