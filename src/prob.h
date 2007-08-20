#include <R.h>
#include <Rmath.h>

#ifndef PROB
#define PROB

typedef struct {
		double **A;
		double **Z;
		double *Pi;
		double *Y;
		double *mu;
		double *sigma;
} params;

params *make_params(int N, int init, int extended);

void params_copy(params *from, params *to, int N);

void print_params(params *p, int N);

double safe_log(double x);

int member(int elem, int array[], int len);

double trans_dist(int dist, double trans, int L, int N);

double deriv_log_trans_dist(int dist, double trans, int L, int N);

double emission_prob(double obs, double mu, double sigma, int p_log);

double prior_prob(params *p, int *Q, int N, int T, int *chrom_starts, int *chroms,
				  int dist, int L, int *distance);

double likelihood_prob(params *p, double *obs, int *Q, int N, int T, int overlap,
					   double *overlaps, int *overlap_ids, int *no_overlaps,
					   int *start_overlaps);

double param_prior_prob(params *p, double *mean_ref, double sd_min, double mean_sd,
						int N, double **W_A, double *W_Pi);

double joint_posterior_prob(params *p, double *obs, int *Q, double *mean_ref,
							int N, int T, double sd_min, double mean_sd,
							int overlap, double *overlaps,
							int *overlap_ids, int *no_overlaps,
							int *start_overlaps, int *chrom_starts,
							int *chroms, int dist, int L, int *distance,
							double **W_A, double *W_Pi);

double Dirichlet(double *x, double *alpha, int N);

double deriv_mu(double mu, double sigma, double obs);

double deriv_sigma(double mu, double sigma, double obs);

double deriv_obs(double mu, double sigma, double ref, double mean_sd);

void printArray(double *x, int len);

void printMatrix(double **x, int rows, int cols);

#endif
