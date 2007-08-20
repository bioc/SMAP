#include "prob.h"

const double H = 0.001;

params *make_params(int N, int init, int extended) {

		params *tmp = (params*) R_alloc(1, sizeof(params));
		
		tmp->A = (double**) R_alloc(N, sizeof(double*));
		tmp->Pi = (double*) R_alloc(N, sizeof(double));
		tmp->mu = (double*) R_alloc(N, sizeof(double));
		tmp->sigma = (double*) R_alloc(N, sizeof(double));

		if (extended) {
				tmp->Z = (double**) R_alloc(N, sizeof(double*));
				tmp->Y = (double*) R_alloc(N, sizeof(double));
		} else {
				tmp->Z = NULL;
				tmp->Y = NULL;
		}
		
		for (int i = 0; i < N; i++) {
				tmp->A[i] = (double*) R_alloc(N, sizeof(double));

				if (extended)
						tmp->Z[i] = (double*) R_alloc(N, sizeof(double));
				
				if (init) {
						for (int j = 0; j < N; j++) {
								tmp->A[i][j] = 0.0;
								if (extended)
										tmp->Z[i][j] = 0.0;
						}
						tmp->Pi[i] = 0.0;
						tmp->mu[i] = 0.0;
						tmp->sigma[i] = 0.0;
						if (extended)
								tmp->Y[i] = 0.0;
				}
		}
		
		return(tmp);
}

void params_copy(params *from, params *to, int N) {

		for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) 
						to->A[i][j] = from->A[i][j];
				to->Pi[i] = from->Pi[i];
				to->mu[i] = from->mu[i];
				to->sigma[i] = from->sigma[i];
		}

		if (from->Z)
				for (int i = 0; i < N; i++) {
						for (int j = 0; j < N; j++)
								to->Z[i][j] = from->Z[i][j];
						if (from->Y)
								to->Y[i] = from->Y[i];
				}
}

void print_params(params *p, int N) {
		
		Rprintf("Pi:\n");
		printArray(p->Pi,N);

		Rprintf("mu:\n");
		printArray(p->mu,N);

		Rprintf("sigma:\n");
		printArray(p->sigma,N);

		Rprintf("A:\n");
		printMatrix(p->A,N,N);	

		if (p->Y) {
				Rprintf("Y:\n");
				printArray(p->Y,N);
		}

		if (p->Z) {
				Rprintf("Z:\n");
				printMatrix(p->Z,N,N);	
		}
}

double safe_log(double x) {

		double x_log = log(x);
		if (!R_finite(x_log)) {
				x_log = -1.0 * DBL_MAX;
				warning("Conversion of log %f to precision in safe_log\n",x);
		}
		
		return(x_log);
}

double trans_dist(int dist, double trans, int L, int N) {

		double td = 0.0;
		if (dist > 0)
				td = trans - ((1.0 - exp(-((double)dist) / ((double)L))) *
							  (trans - (((double)1) / ((double)N))));
		else
				td = trans;

		return(td);
}

double deriv_log_trans_dist(int dist, double trans, int L, int N) {
		
		if (dist > 0)
				return((1.0 / trans_dist(dist, trans, L, N)) *
					   (exp(-((double)dist) / ((double)L))));
		else
				return(1.0 / trans);
}

int member(int elem, int array[], int len) {
		for (int i = 0; i < len; i++) {
				if (elem == array[i])
						return(1);
		}
		return(0);
}

double emission_prob(double obs, double mu, double sigma, int p_log) {
		
		int lower = obs < mu;
		double p = pnorm(obs, mu, sigma, lower, p_log);
		if (p_log && !R_finite(p)) {
				p = -1.0 * DBL_MAX;
				warning("Conversion of %f to precision in emission_prob\n",p);
		}
		if (p_log)
				p += M_LN2;
		else
				p = 2 * p;
		return(p);
}

double prior_prob(params *p, int *Q, int N, int T, int *chrom_starts, int *chroms,
				  int dist, int L, int *distance) {

		double prob = 0.0;

		for (int i = 0; i < *chroms; i++) {

				prob += safe_log(p->Pi[Q[chrom_starts[i]]]);

				int start = chrom_starts[i];
				int end = T-1;
				if (i < (*chroms - 1))
						end = chrom_starts[i+1]-1;
				
				if ((end - start) > 0) {
						for (int t = start; t < end; t++) {
								if (dist)
										prob += safe_log(trans_dist(distance[t+1],
																	p->A[Q[t]][Q[t+1]],
																	L, N));
								else
										prob += safe_log(p->A[Q[t]][Q[t+1]]);
						}
				}
		}

		return(prob);
}

double likelihood_prob(params *p, double *obs, int *Q, int N, int T, int overlap,
					   double *overlaps, int *overlap_ids, int *no_overlaps,
					   int *start_overlaps) {

		double prob = emission_prob(obs[0], p->mu[Q[0]], p->sigma[Q[0]], 1);

		if (T > 1) {
				int no_olaps;
				int start;
				double sum_olap;
				for (int t = 1; t < T; t++) {
						if (overlap) {
								no_olaps = no_overlaps[t];
								int olap_ids[no_olaps];
								double olaps[no_olaps];
								int qts[no_olaps];
								start = start_overlaps[t];
								sum_olap = 1.0;
								for (int i = 0; i < no_olaps; i++) {
										olap_ids[i] = overlap_ids[start + i];
										qts[i] = Q[olap_ids[i]];
										olaps[i] = overlaps[start + i];
										sum_olap += olaps[i];
								}
								olaps[no_olaps-1] = 1.0;							
								

								for (int i = 0; i < no_olaps; i++) {
										prob += emission_prob(obs[t], p->mu[qts[i]],
															  p->sigma[qts[i]], 1) +
												safe_log(olaps[i] / sum_olap);
								}
						}
						else {
								prob += emission_prob(obs[t], p->mu[Q[t]], p->sigma[Q[t]], 1);
						}
				}
		}
		return(prob);
}

double param_prior_prob(params *p, double *mean_ref, double sd_min, double mean_sd,
						int N, double **W_A, double *W_Pi) {

		double prob = 0.0;

		for (int i = 0; i < N; i++) {
				prob += safe_log(Dirichlet(p->A[i], W_A[i], N));
				prob += safe_log(sd_min / p->sigma[i]) +
						emission_prob(p->mu[i], mean_ref[i], mean_sd, 1);
		}
		prob += safe_log(Dirichlet(p->Pi, W_Pi, N));

		return(prob);
}

double Dirichlet(double *x, double *alpha, int N) {

		double p = 1.0;
		double numer = 1.0;
		double denom = 0.0;

		for (int i = 0; i < N; i++) {

				numer *= gammafn(alpha[i]);
				denom += alpha[i];
				p *= R_pow(x[i], alpha[i] - 1);
		}
		
		denom = gammafn(denom);
		p *= numer / denom;
		
		return(p);
}

double joint_posterior_prob(params *p, double *obs, int *Q, double *mean_ref,
							int N, int T, double sd_min, double mean_sd,
							int overlap, double *overlaps,
							int *overlap_ids, int *no_overlaps,
							int *start_overlaps, int *chrom_starts, int *chroms,
							int dist, int L, int *distance, double **W_A, double *W_Pi) {

		double prob = prior_prob(p, Q, N, T, chrom_starts, chroms, dist, L, distance);
		prob += likelihood_prob(p, obs, Q, N, T, overlap, overlaps, overlap_ids,
								no_overlaps, start_overlaps);
		prob += param_prior_prob(p, mean_ref, sd_min, mean_sd, N, W_A, W_Pi);
		
		return(prob);
}

double deriv_mu(double mu, double sigma, double obs) {

		return((emission_prob(obs, mu+H, sigma, 0) -
				emission_prob(obs, mu, sigma, 0)) / H);
}

double deriv_sigma(double mu, double sigma, double obs) {

		return((emission_prob(obs, mu, sigma+H, 0) -
				emission_prob(obs, mu, sigma, 0)) / H);
}

double deriv_obs(double mu, double sigma, double ref, double mean_sd) {

		return(((emission_prob(mu+H, ref, sigma, 0) -
				 emission_prob(mu, ref, sigma, 0)) / H) /
			   emission_prob(mu, ref, mean_sd, 0));
}


void printArray(double *x, int len) {

		for (int i = 0; i < len; i++) {
				Rprintf("%f ", x[i]);
		}
		Rprintf("\n");
}

void printMatrix(double **x, int rows, int cols) {
		
		for (int i = 0; i < rows; i++) {
				Rprintf("%d: ", i);
				printArray(x[i], cols);
		}
}

/*
void printIntArray(int *x, int len) {

		for (int i = 0; i < len; i++) {
				Rprintf("%d ", x[i]);
		}
		Rprintf("\n");
}

void printIntMatrix(int *x, int rows, int cols) {

		for (int i = 0; i < rows; i++) {
				Rprintf("%d: ", i);
				printIntArray(&x[i], cols);
		}
}
*/
