#include "overlap.h"

o_struct *make_o_struct(int n) {

		o_struct *tmp = (o_struct*) R_alloc(1, sizeof(o_struct));

		tmp->n = n;
		tmp->start = 0;
		tmp->ids = (int*) R_alloc(n, sizeof(int));
		tmp->olaps = (double*) R_alloc(n, sizeof(double));

		return(tmp);
}

void print_struct(o_struct *o, int n) {

		Rprintf("n: %d\n", o->n);
		Rprintf("start: %d\n", o->start);
		Rprintf("ids: ");
		for (int i = 0; i < n; i++)
				Rprintf("%d ", o->ids[i]);
		Rprintf("\n");
		Rprintf("olaps: ");
		for (int i = 0; i < n; i++)
				Rprintf("%f ", o->olaps[i]);
		Rprintf("\n");
}

double calc_overlap(double s1, double e1, int c1, double s2, double e2, int c2) {

		if ((c1 != c2) || (e2 == s2))
				return(0.0);

		return(min( max( (min(e1, e2) - max(s1, s2)) / (e2 - s2), 0.0), 1.0));
}

SEXP calc_overlaps(SEXP T, SEXP startpos, SEXP endpos, SEXP chrom) {

		SEXP res, olap_ids, olaps, start_olaps, no_olaps;

		int t = INTEGER(T)[0];
		
		double *d_startpos = REAL(startpos);
		double *d_endpos = REAL(endpos);
		int *i_chrom = INTEGER(chrom);
		
		o_struct **overlap = (o_struct**) R_alloc(t, sizeof(o_struct*));

		overlap[0] = make_o_struct(1);
		overlap[0]->start = 1;
		overlap[0]->ids[0] = 1;
		overlap[0]->olaps[0] = 0.0;

		int prev_zero = 0;
		int start = 2;
		int no_ids;
		int iter;

		int len = 1;
		
		if (t > 1) {
				for (int i = 1; i < t; i++) {
						
						double *olap = (double*) R_alloc(i - prev_zero, sizeof(double));
						iter = 0;
						no_ids = 0;
						for (int j = prev_zero; j < i; j++) {
								olap[iter] = calc_overlap(d_startpos[j], d_endpos[j], i_chrom[j],
														  d_startpos[i], d_endpos[i], i_chrom[i]);
								if (olap[iter] > 0.0)
										no_ids++;
								iter++;
						}
						
						no_ids++;
						overlap[i] = make_o_struct(no_ids);
						overlap[i]->start = start;
						start = start + no_ids;

						len += no_ids;

						if (no_ids > 1) {
								iter = 0;
								for (int j = 0; j < (i - prev_zero); j++) {
										if (olap[j] > 0.0) {
												overlap[i]->ids[iter] = j + prev_zero + 1;
												overlap[i]->olaps[iter] = olap[j];
												iter++;
										}
								}
						}
						else {
								prev_zero = i;
						}
						
						overlap[i]->ids[no_ids-1] = i+1;
						overlap[i]->olaps[no_ids-1] = 0.0;
				}
		}
		
		PROTECT(start_olaps = allocVector(INTSXP, t));
		PROTECT(no_olaps = allocVector(INTSXP, t));
		PROTECT(olaps = allocVector(REALSXP, len));
		PROTECT(olap_ids = allocVector(INTSXP, len));

		int *i_start_olaps = INTEGER(start_olaps);
		int *i_no_olaps = INTEGER(no_olaps);
		double *d_olaps = REAL(olaps);
		int *i_olap_ids = INTEGER(olap_ids);

		int len_iter = 0;
		o_struct *o;
		for (int i = 0; i < t; i ++) {

				o = overlap[i];
				i_start_olaps[i] = o->start;
				i_no_olaps[i] = o->n;

				for (int j = 0; j < o->n; j++, len_iter++) {
						d_olaps[len_iter] = o->olaps[j];
						i_olap_ids[len_iter] = o->ids[j];
				}
		}
		
		PROTECT(res = allocVector(VECSXP, 4));
		SET_VECTOR_ELT(res, 0, olap_ids);
		SET_VECTOR_ELT(res, 1, olaps);
		SET_VECTOR_ELT(res, 2, start_olaps);
		SET_VECTOR_ELT(res, 3, no_olaps);
		
		UNPROTECT(5);
		
		return(res);
}
