/* includes routines to solve linear equations systems, given that the matrix of the coefficients is symmetric.
	the functions can input the entries of the symmetric matrix as a bidimensional as well as an unidimensional array. */

#ifndef LLT_ROUTINES
#define LLT_ROUTINES

#include <math.h>


using namespace std;


/** Given a matrix a[0..n-1][0..n-1], this routine replaces it by its L*Lt decomposition.
a and n are input. a is output, arranged as follows:

	[ @00	]
	[ @10 @11	]
	[ @20 @21 @22	]
	[ @30 @31 @32 @33 ]

This routine is used together with <code>LLt_backsubst(double *(&a), long long int n, double *(&b))</code>
to solve linear equations systems. */

void LLt_factorize (double *(&a), long int n);

/** Does the same as <code>LLt_factorize(double *(&a), long long int n)</code> only the matrix is a bidimensional array.
This routine is used together with <code>LLt_backsubst (double **(&a), long long int n, double *(&b))</code>
to solve linear equations systems. */

void LLt_factorize (double **(&a), long int n);


/*
Solves the set of n linear equations A * x = b. Here a[0..n-1][0..n-1] is input, not as the matrix
A but rather as its L*Lt decompostion, determined by the routine LLt_factorize. b[0..n-1] is input
as the right-hand side vector b, and returns with the solution vector x. a and n are not modified
by this routine.
*/

void LLt_backsubst (double *(&a), long int n, double *(&b));

void LLt_backsubst (double **(&a), long int n, double *(&b));


#endif



