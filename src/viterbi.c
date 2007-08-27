#include "viterbi.h"

void viterbi(int *x_T, int *x_N, double *x_A, double *x_Pi, double *mu, double *sigma,
			 double *obs, int *overlap, double *overlaps, int *overlap_ids, int *no_overlaps,
			 int *start_overlaps, int *dist, int *L, int *distance, double *P, int *Q,
			 double *mean_ref, double *sd_min, double *mean_sd, int *prior,
			 double *x_W_A, double *W_Pi) {
		
		int N = *x_N;
		int T = *x_T;
		double A[N][N];
		double W_A[N][N];
		double Pi[N];
		
		double delta[T][N];
		int psi[T][N];
		
		// Fill A and Pi
		for (int i = 0; i < N; i++) {
				for (int j = 0, index = i; j < N; j++, index += N) {
						if (*dist)
								A[i][j] = x_A[index];
						else
								A[i][j] = safe_log(x_A[index]);
						W_A[i][j] = x_W_A[index];
				}
				Pi[i] = safe_log(x_Pi[i]);
		}
				
		// Initialization
		for (int i = 0; i < N; i++) {
				delta[0][i] = emission_prob(obs[0], mu[i], sigma[i], 1) + Pi[i];
		}
		
		// Recursion
		int no_olaps;
		int start;
		double sum_olap;
		double trans;
		if (T > 1) {
				for (int t = 1; t < T; t++) {
						
						no_olaps = no_overlaps[t];
						int olap_ids[no_olaps];
						double olaps[no_olaps];
						start = start_overlaps[t];
						sum_olap = 1.0;
						
						if (*overlap) {
								for (int i = 0; i < no_olaps; i++) {
										olap_ids[i] = overlap_ids[start + i];
										olaps[i] = overlaps[start + i];
										sum_olap += olaps[i];
								}
								olaps[no_olaps-1] = 1.0;
						}
						
						for (int j = 0; j < N; j++) {
								double prev[N];
								if (*dist)
										prev[0] = safe_log(trans_dist(distance[t], A[0][j],
																	  *L, N)) +
												delta[t-1][0];
								else
										prev[0] = A[0][j] + delta[t-1][0];
								double max = prev[0];
								int maxid = 0;
								if (N > 1) {
										for (int i = 1; i < N; i++) {
												if (*dist)
														prev[i] = safe_log(trans_dist(distance[t],
																					  A[i][j],
																					  *L, N)) +
																delta[t-1][i];
												else
														prev[i] = A[i][j] + delta[t-1][i];
												if (prev[i] > max) {
														maxid = i;
														max = prev[i];
												}
										}
								}
								
								psi[t][j] = maxid;

								trans = 0.0;
								
								if (*overlap) {
										int qt[no_olaps];
										if (no_olaps > 1) {
												int q = j;
												int iter = no_olaps-2;
												for (int i = t-1; i >= olap_ids[0]; i--) {
														q = psi[i+1][q];
														if (member(i, olap_ids, no_olaps)) {
																qt[iter] = q;
																iter--;
														}
												}
										}
										qt[no_olaps-1] = j;
										
										int id;
										for (int i = 0; i < no_olaps; i++) {
												id = qt[i];
												trans += emission_prob(obs[t], mu[id], sigma[id], 1) +
														safe_log(olaps[i] / sum_olap);
										}
								}
								else {
										trans = emission_prob(obs[t], mu[j], sigma[j], 1);
								}

								if (*dist)
										delta[t][j] = delta[t-1][psi[t][j]] +
												safe_log(trans_dist(distance[t], A[psi[t][j]][j], *L, N)) +
												trans;
								else
										delta[t][j] = delta[t-1][psi[t][j]] + A[psi[t][j]][j] + trans;
						}
				}
		}
		
		// Termination
		double max = delta[T-1][0];
		int maxid = 0;
		if (N > 1) {
				for (int i = 1; i < N; i++) {
						if (delta[T-1][i] > max) {
								maxid = i;
								max = delta[T-1][i];
						}
				}
		}
		Q[T-1] = maxid;
		*P = delta[T-1][Q[T-1]];

		// Calculate parameter prior probability
		if (*prior) {
				for (int i = 0; i < N; i++) {
						*P += safe_log(Dirichlet(A[i], W_A[i], N));
						*P += safe_log(*sd_min / sigma[i]) +
								emission_prob(mu[i], mean_ref[i], *mean_sd, 1);
				}
				*P += safe_log(Dirichlet(Pi, W_Pi, N));
		}
		
		// Path backtracking
		if (T > 1) {
				for (int t = T-2; t >= 0; t--) {
						Q[t] = psi[t+1][Q[t+1]];
				}
		}
}
