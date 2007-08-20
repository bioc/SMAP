#include <R.h>
#include <Rinternals.h>
#include "overlap.h"
#include "gradient.h"
#include "viterbi.h"
#include <R_ext/Rdynload.h>

static R_NativePrimitiveArgType gradient_descent_t[37] = {
		INTSXP,  /* 		int *_T, */
		INTSXP,  /* 		int *_N, */
		REALSXP, /* 		double *Pi, */
		REALSXP, /* 		double *_A, */
		REALSXP, /* 		double *Y, */
		REALSXP, /* 		double *_Z, */
		REALSXP, /* 		double *mu, */
		REALSXP, /* 		double *sigma, */
		INTSXP,  /* 		int *Q, */
		REALSXP, /* 		double *_P, */
		INTSXP,  /* 		int *len, */
		REALSXP, /* 		double *obs, */
		INTSXP,  /* 		int *overlap, */
		REALSXP, /* 		double *overlaps, */
		INTSXP,  /* 		int *overlap_ids, */
		INTSXP,  /* 		int *no_overlaps, */
		INTSXP,  /* 		int *start_overlaps, */
		INTSXP,  /* 		int *dist, */
		INTSXP,  /* 		int *L, */
		INTSXP,  /* 		int *distance, */
		INTSXP,  /* 		int *chrom_starts, */
		INTSXP,  /* 		int *chroms, */
		REALSXP, /* 		double *mean_ref, */
		REALSXP, /* 		double *sd_min, */
		REALSXP, /* 		double *mean_sd, */
		INTSXP,  /* 		int *max_iters, */
		INTSXP,  /* 		int *inf_iters, */
		REALSXP, /* 		double *_tau, */
		REALSXP, /* 		double *eta, */
		REALSXP, /* 		double *e_change, */
		REALSXP, /* 		double *e_same, */
		REALSXP, /* 		double *e_min, */
		REALSXP, /* 		double *e_max, */
		INTSXP,  /* 		int *adaptive, */
		INTSXP,  /* 		int *verbose */
		REALSXP, /* 		double *W_A, */
		REALSXP  /* 		double *W_Pi */
};

static R_NativePrimitiveArgType viterbi_t[23] = {
		INTSXP,  /* 		int *_T, */
		INTSXP,  /* 		int *_N, */
		REALSXP, /* 		double *_A, */
		REALSXP, /* 		double *_Pi, */
		REALSXP, /* 		double *mu, */
		REALSXP, /* 		double *sigma, */
		REALSXP, /* 		double *obs, */
		INTSXP,  /* 		int *overlap, */
		REALSXP, /* 		double *overlaps, */
		INTSXP,  /* 		int *overlap_ids, */
		INTSXP,  /* 		int *no_overlaps, */
		INTSXP,  /* 		int *start_overlaps, */
		INTSXP,  /* 		int *dist, */
		INTSXP,  /* 		int *L, */
		INTSXP,  /* 		int *distance, */
		REALSXP, /* 		double *P, */
		INTSXP,  /* 		int *Q, */
		REALSXP, /* 		double *mean_ref, */
		REALSXP, /* 		double *sd_min, */
		REALSXP, /* 		double *mean_sd, */
		INTSXP,  /* 		int *prior */
		REALSXP, /* 		double *W_A, */
		REALSXP  /* 		double *W_Pi */
};

static const R_CMethodDef CEntries[]  = {
		{"gradient_descent", (DL_FUNC) &gradient_descent, 37, gradient_descent_t},
		{"viterbi", (DL_FUNC) &viterbi, 23, viterbi_t},
		{NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
		{"calc_overlaps", (DL_FUNC) &calc_overlaps, 4},
		{NULL, NULL, 0}
};

void
R_init_SMAP(DllInfo *dll)
{
		R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
		R_useDynamicSymbols(dll, FALSE);
}
