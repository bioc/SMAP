#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include "prob.h"

#ifndef GRADIENT
#define GRADIENT

void gradient(int start, int T, int N, params *grad, params *hmm,
			  int *Q, int *len, double *obs, int *overlap,
			  double *overlaps, int *overlap_ids, int *no_overlaps,
			  int *start_overlaps, int dist, int L, int *distance);

void prior_gradient(params *grad, params *hmm, int N, double *mean_ref,
					double *sd_min, double *mean_sd, double **W_A, double *W_Pi);

void normalize(params *grad, int N);

void scale_eta(params *eta, double e_change, int N);

void eta_change(double *eta, double *grad, int sign, double e_change, double e_same,
				double e_min, double e_max);

void eta_update(params *eta, params *grad, params *work_grad,
				double e_change, double e_same, double e_min, double e_max, int N);

void hmm_update(params *hmm, params *eta, params *grad, int N, double sd_min);

void gradient_descent(int *x_T, int *x_N, double *Pi, double *x_A,
					  double *Y, double *x_Z,
					  double *mu, double *sigma, int *Q,
					  double *x_P, int *len, double *obs, int *overlap, double *overlaps,
					  int *overlap_ids, int *no_overlaps,
					  int *start_overlaps, int *dist, int *L, int *distance, int *chrom_starts,
					  int *chroms, double *mean_ref,
					  double *sd_min, double *mean_sd,
					  int *max_iters, int *inf_iters, double *x_tau,
					  double *eta, double *e_change,
					  double *e_same, double *e_min, double *e_max,
					  int *adaptive, int *verbose, double *x_W_A, double *W_Pi);

#endif
