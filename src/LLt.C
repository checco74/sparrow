
#include "LLt.h"


using namespace std;




void LLt_factorize (double *(&a), long int n)
{
	long int i, j, k, ij, ic = 0, jc = 0;
	long int ik, jk, jj;
	double sum = 0.0;
	for (j = 0; j<n; j++)
	{
	jc += j; ic = jc;
	for (k = 0; k<j; k++) {jk = jc+k; sum += a[jk]*a[jk];}
	jj = jc+j; a[jj] = sqrt(fabs(a[jj]-sum));
	sum = 0.0;
	for (i = j+1; i<n; i++)
	{
	ic += i;
	for (k = 0; k<j; k++) {ik = ic+k; jk = jc+k; sum += a[ik]*a[jk];}
	ij = ic+j; a[ij] = (a[ij]-sum)/a[jj];
	sum = 0.0;
	}
	}
}

void LLt_factorize (double **(&a), long int n)
{
	long int i, j, k;
	double sum = 0.0;
	for (j = 0; j<n; j++)
	{
	for (k = 0; k<j; k++) {sum += a[j][k]*a[j][k];}
	a[j][j] = sqrt(fabs(a[j][j]-sum));
	sum = 0.0;
	for (i = j+1; i<n; i++)
	{
	for (k = 0; k<j; k++) {sum += a[i][k]*a[j][k];}
	a[i][j] = (a[i][j]-sum)/a[j][j];
	sum = 0.0;
	}
	}
}







void LLt_backsubst (double *(&a), long int n, double *(&b))
{
	double sum = 0.0;
	long int i, j, ic = 0, ij = 0;
	long int nn = (n*(n-1))/2;
	for (i = 0; i<n; i++)
	{
	ic += i;
	sum = 0.0;
	for (j = 0; j<i; j++) {ij = ic+j; sum += a[ij]*b[j];}
	b[i] = (b[i]-sum)/a[ic+i];
	}
	ic = nn;
	for (i = n-1; i>=0; i--)
	{
	ij = i+nn;
	sum = 0.0;
	for (j = n-1; j>i; j--) {sum += a[ij]*b[j]; ij -= j;}
	b[i] = (b[i]-sum)/a[ic+i];
	ic -= i;
	}
}

void LLt_backsubst (double **(&a), long int n, double *(&b))
{
	long int i, j;
	double sum = 0.0;
	for (i = 0; i<n; i++)
	{
	sum = 0.0;
	for (j = 0; j<i; j++) sum += a[i][j]*b[j];
	b[i] = (b[i]-sum)/a[i][i];
	}
	for (i = n-1; i>=0; i--)
	{
	sum = 0.0;
	for (j = n-1; j>i; j--) sum += a[j][i]*b[j];
	b[i] = (b[i]-sum)/a[i][i];
	}
}




