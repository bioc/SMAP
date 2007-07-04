#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

#ifndef OVERLAP
#define OVERLAP

#define min( A , B )  ( ( A ) < ( B ) ? ( A ) : ( B ) )
#define max( A , B )  ( ( A ) > ( B ) ? ( A ) : ( B ) )

typedef struct {
		int n;
		int start;
		int *ids;
		double *olaps;
} o_struct;

o_struct *make_o_struct(int n);

void print_struct(o_struct *o, int n);

double calc_overlap(double s1, double e1, int c1, double s2, double e2, int c2);

SEXP calc_overlaps(SEXP T, SEXP startpos, SEXP endpos, SEXP chrom);

#endif
