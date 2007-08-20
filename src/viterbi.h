#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include "prob.h"

#ifndef VITERBI
#define VITERBI

void viterbi(int *x_T, int *x_N, double *x_A, double *x_Pi, double *mu, double *sigma,
			 double *obs, int *overlap, double *overlaps, int *overlap_ids, int *no_overlaps,
			 int *start_overlaps, int *dist, int *L, int *distance, double *P, int *Q,
			 double *mean_ref, double *sd_min, double *mean_sd, int *prior,
			 double *x_W_A, double *W_Pi);

#endif
