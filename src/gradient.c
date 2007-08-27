#include "gradient.h"

void gradient(int start_index, int T, int N, params *grad, params *hmm,
			  int *Q, int *len, double *obs, int *overlap,
			  double *overlaps, int *overlap_ids, int *no_overlaps,
			  int *start_overlaps, int dist, int L, int *distance) {
		
		// Pi
		grad->Pi[Q[start_index]] -= 1 / hmm->Y[Q[start_index]];
		
		// A
		if ((T - start_index) > 0) {
				int i, j;
				for (int t = start_index; t < T; t++) {
						i = Q[t];
						j = Q[t+1];
						if (dist)
								grad->A[i][j] -= deriv_log_trans_dist(distance[t+1], hmm->Z[i][j],
																	  L, N);
						else
								grad->A[i][j] -= 1 / hmm->Z[i][j];
				}
		}

		// Phi (obs 1)
		int qt = Q[start_index];
		
		double p = emission_prob(obs[start_index], hmm->mu[qt], hmm->sigma[qt], 0);
		
		// mu
		double mu_deriv = deriv_mu(hmm->mu[qt], hmm->sigma[qt], obs[start_index]) / p;
		
		if (!ISNAN(mu_deriv))
				grad->mu[qt] -= mu_deriv;
		//else
		//warning("mu_deriv is NAN\n");

		//sigma
		double weight = 0.0;
		double sigma_deriv;
		if (len[qt] != 0) {
				weight = 1.0 / (double)len[qt];
				
				sigma_deriv = weight * deriv_sigma(hmm->mu[qt], hmm->sigma[qt], obs[start_index]) / p;

				if (!ISNAN(sigma_deriv))
						grad->sigma[qt] -= sigma_deriv;
				//else
				//warning("sigma_deriv is NAN\n");
		}

		// Phi (obs t)
		int no_olaps;
		int start;
		double sum_olap;
		double trans;

		if ((T - start_index) > 0) {
				for (int t = start_index+1; t <= T; t++) {

						if (*overlap) {
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
								
								trans = 0.0;
								
								for (int i = 0; i < no_olaps; i++) {
										trans += emission_prob(obs[t], hmm->mu[qts[i]],
															   hmm->sigma[qts[i]], 0) *
												olaps[i] / sum_olap;
								}
								
								// mu				
								for (int i = 0; i < no_olaps; i++) {
										mu_deriv = deriv_mu(hmm->mu[qts[i]], hmm->sigma[qts[i]],
															obs[t]) *
												(olaps[i] / sum_olap) / trans;
										if (!ISNAN(mu_deriv))
												grad->mu[qts[i]] -= mu_deriv;
										//else
										//warning("mu_deriv is NAN\n");
								}
								
								// sigma
								for (int i = 0; i < no_olaps; i++) {
										if (len[qts[i]] != 0) {
												weight = 1.0 / (double)len[qts[i]];
												sigma_deriv = weight *
														deriv_sigma(hmm->mu[qts[i]],
																	hmm->sigma[qts[i]], obs[t]) *
														(olaps[i] / sum_olap) / trans;
												if (!ISNAN(sigma_deriv))
														grad->sigma[qts[i]] -= sigma_deriv;
												//else
												//warning("sigma_deriv is NAN\n");
										}
								}
						}
						else {
								qt = Q[t];

								//mu
								p = emission_prob(obs[t], hmm->mu[qt], hmm->sigma[qt], 0);
								mu_deriv = deriv_mu(hmm->mu[qt], hmm->sigma[qt], obs[t]) / p;

								if (!ISNAN(mu_deriv))
										grad->mu[qt] -= mu_deriv;
								//else
								//warning("mu_deriv is NAN\n");
								
								//sigma
								if (len[qt] != 0) {
										weight = 1.0 / (double)len[qt];
										sigma_deriv = weight *
												deriv_sigma(hmm->mu[qt],
															hmm->sigma[qt], obs[t]) / p;

										if (!ISNAN(sigma_deriv))
												grad->sigma[qt] -= sigma_deriv;
										//else 
										//warning("sigma_deriv is NAN\n");
								}
						}
				}
		}		
}

void prior_gradient(params *grad, params *hmm, int N, double *mean_ref,
					double *sd_min, double *mean_sd, double **W_A, double *W_Pi) {

		int lower;
		double mu_deriv_obs;
		
		for (int i = 0; i < N; i++) {

				//A
				for (int j = 0; j < N; j++)
						grad->A[i][j] -= (W_A[i][j] - 1) / hmm->Z[i][j];

				//Pi
				grad->Pi[i] -= (W_Pi[i] - 1) / hmm->Y[i];
				
				lower = hmm->mu[i] < mean_ref[i];
				
				//mu
				mu_deriv_obs = deriv_obs(hmm->mu[i], hmm->sigma[i], mean_ref[i], *mean_sd);
				if (!ISNAN(mu_deriv_obs))
						grad->mu[i] -= mu_deriv_obs;
				//else
				//warning("mu_deriv_obs is NAN\n");
				
				//sigma
				if (hmm->sigma[i] > *sd_min)
						grad->sigma[i] += 1 / hmm->sigma[i];
		}
}

void normalize(params *grad, int N) {
		double norm = 0.0;
		for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++)
						norm += R_pow_di(grad->A[i][j], 2);
				norm += R_pow_di(grad->Pi[i], 2);
				norm += R_pow_di(grad->mu[i], 2);
				norm += R_pow_di(grad->sigma[i], 2);
		}
		norm = sqrt(norm);
		for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++)
						grad->A[i][j] /= norm;
				grad->Pi[i] /= norm;
				grad->mu[i] /= norm;
				grad->sigma[i] /= norm;
		}
}

void scale_eta(params *eta, double e_change, int N) {

		for (int i = 0; i < N; i++) {
				//Update Pi
				eta->Pi[i] *= e_change;
				// Update A
				for (int j = 0; j < N; j++)
						eta->A[i][j] *= e_change;
				// Update mu
				eta->mu[i] *= e_change;
				// Update sigma
				eta->sigma[i] *= e_change;
		}
}

void eta_change(double *eta, double *grad, int sign, double e_change, double e_same,
				double e_min, double e_max) {

		if (sign == -1) {
				(*eta) *= e_change;
				if ((*eta) < e_min)
						(*eta) = e_min;
				(*grad) = 0.0;
		}
		else if (sign == 1) {
				(*eta) *= e_same;
				if ((*eta) > e_max)
						(*eta) = e_max;
		}
}

void eta_update(params *eta, params *grad, params *work_grad,
				double e_change, double e_same, double e_min, double e_max, int N) {

		double prod;
		
		for (int i = 0; i < N; i++) {
				// Update Pi
				prod = grad->Pi[i] * work_grad->Pi[i];
				eta_change(&eta->Pi[i], &work_grad->Pi[i], sign(prod), e_change, e_same, e_min,
						   e_max);
				// Update A
				for (int j = 0; j < N; j++) {
						prod = grad->A[i][j] * work_grad->A[i][j];
						eta_change(&eta->A[i][j], &work_grad->A[i][j], sign(prod), e_change, e_same,
								   e_min, e_max);
				}
				// Update mu
				prod = grad->mu[i] * work_grad->mu[i];
				eta_change(&eta->mu[i], &work_grad->mu[i], sign(prod), e_change, e_same, e_min, e_max);
				// Update sigma
				prod = grad->sigma[i] * work_grad->sigma[i];
				eta_change(&eta->sigma[i], &work_grad->sigma[i], sign(prod), e_change, e_same,
						   e_min, e_max);
		}
}

void hmm_update(params *hmm, params *eta, params *grad, int N, double sd_min) {

		double tmp;
		double Pi_sum = 0.0;
		double A_sum[N];
		for (int i = 0; i < N; i++) {
				// Update Pi and Y
				hmm->Pi[i] *= exp(-1.0 * eta->Pi[i] * grad->Pi[i]);
				hmm->Y[i] -= eta->Pi[i] * grad->Pi[i];
				Pi_sum += hmm->Pi[i];

				// Update A and Z
				A_sum[i] = 0.0;
				for (int j = 0; j < N; j++) {
						hmm->A[i][j] *= exp(-1.0 * eta->A[i][j] * grad->A[i][j]);
						hmm->Z[i][j] -= eta->A[i][j] * grad->A[i][j];
						A_sum[i] += hmm->A[i][j];
				}

				// Update mu
				hmm->mu[i] -= eta->mu[i] * grad->mu[i];
				
				// Update sigma
				hmm->sigma[i] -= eta->sigma[i] * grad->sigma[i];
				if (hmm->sigma[i] <= sd_min) {
						hmm->sigma[i] = sd_min;
						grad->sigma[i] = 0.0;
				}
		}

		//Normalize Pi and A
		for (int i = 0; i < N; i++) {
				hmm->Pi[i] /= Pi_sum;
				for (int j = 0; j < N; j++)
						hmm->A[i][j] /= A_sum[i];
		}
}

void gradient_descent(int *x_T, int *x_N, double *Pi, double *x_A,
					  double *Y, double *x_Z,
					  double *mu, double *sigma, int *Q,
					  double *x_P, int *len, double *obs, int *overlap, double *overlaps,
					  int *overlap_ids, int *no_overlaps,
					  int *start_overlaps, int *dist, int *L, int *distance,
					  int *chrom_starts,
					  int *chroms, double *mean_ref,
					  double *sd_min, double *mean_sd,
					  int *max_iters, int *inf_iters, double *x_tau,
					  double *eta, double *e_change,
					  double *e_same, double *e_min, double *e_max,
					  int *adaptive, int *verbose, double *x_W_A, double *W_Pi) {
		
		int N = *x_N;
		int T = *x_T;
		double P = *x_P;
		double tau = *x_tau;
		
		double **W_A = (double**) R_alloc(N, sizeof(double*));
		
		params *grad_params = make_params(N,1,0);
		params *work_grad_params = make_params(N,1,0);
		params *eta_params = make_params(N,0,0);
		params *hmm_params = make_params(N,0,1);
		params *work_hmm_params = make_params(N,0,1);

		// Initialize A, Z, hmm_params
		for (int i = 0; i < N; i++) {
				W_A[i] = (double*) R_alloc(N, sizeof(double));
				for (int j = 0, index = i; j < N; j++, index += N) {
						hmm_params->A[i][j] = x_A[index];
						hmm_params->Z[i][j] = x_Z[index];
						W_A[i][j] = x_W_A[index];
						eta_params->A[i][j] = *eta;
				}
				hmm_params->Y[i] = Y[i];
				hmm_params->Pi[i] = Pi[i];
				hmm_params->mu[i] = mu[i];
				hmm_params->sigma[i] = sigma[i];
				eta_params->Pi[i] = *eta;
				eta_params->mu[i] = *eta;
				eta_params->sigma[i] = *eta;
		}

		params_copy(hmm_params, work_hmm_params, N);
		
		int opt = 0;
		int iter = 1;
		int zero_grad = 1;
			
		double new_P = P;

		while (!opt && (*inf_iters || (iter <= *max_iters))) {

				R_CheckUserInterrupt();
				
				// Calculate gradient
				for (int i = 0; i < N; i++) {
						for (int j = 0; j < N; j++) {
								work_grad_params->A[i][j] = 0.0;
						}
						work_grad_params->Pi[i] = 0.0;
						work_grad_params->mu[i] = 0.0;
						work_grad_params->sigma[i] = 0.0;
				}

				int start;
				int end;
				for (int i = 0; i < *chroms; i++) {
						start = chrom_starts[i];
						end = T-1;
						if (i < (*chroms - 1))
								end = chrom_starts[i+1]-1;
						
						gradient(start, end, N, work_grad_params, work_hmm_params, Q, len, obs,
								 overlap, overlaps, overlap_ids, no_overlaps,
								 start_overlaps, *dist, *L, distance);
				}

				// Calculate gradient for prior
				prior_gradient(work_grad_params, work_hmm_params, N, mean_ref, sd_min, mean_sd,
							   W_A, W_Pi);
				
				// Normalize gradient
				normalize(work_grad_params, N);

				// Update eta from orig grad and new grad
				if (*adaptive && !zero_grad)
						eta_update(eta_params, grad_params, work_grad_params,
								   *e_change, *e_same, *e_min, *e_max, N);
				
				params_copy(work_grad_params, grad_params, N);
				
				// Update working hmm parameters
				hmm_update(work_hmm_params, eta_params, grad_params, N, *sd_min);
				
				// Calculate joint posterior prob
				new_P = joint_posterior_prob(work_hmm_params, obs, Q, mean_ref, N, T,
											 *sd_min, *mean_sd, *overlap,
											 overlaps, overlap_ids, no_overlaps, start_overlaps,
											 chrom_starts, chroms, *dist, *L, distance,
											 W_A, W_Pi);
				
				// Improvement?
				if (new_P > P + tau) {
						if (*verbose > 2) {
								Rprintf("   GD iteration %d, P: %f\n", iter, new_P);
						}
						
						// Update hmm from working params
						params_copy(work_hmm_params, hmm_params, N);
						
						P = new_P;
						zero_grad = 0;
				} else {
						opt = 1;
				}
				
				iter++;
		}

		*x_P = P;
		
		//Update Pi, x_A, x_Z, Y, mu, sigma
		int iterator = 0;
		for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
						x_A[iterator] = hmm_params->A[j][i];
						x_Z[iterator] = hmm_params->Z[j][i];
						iterator++;
				}
				Pi[i] = hmm_params->Pi[i];
				Y[i] = hmm_params->Y[i];
				mu[i] = hmm_params->mu[i];
				sigma[i] = hmm_params->sigma[i];
		}
}
