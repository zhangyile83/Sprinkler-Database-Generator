#include <conio.h>
#include <stdlib.h>
#include <math.h>

//polyfit(n,x,y,poly_n,a)
//���y=a0+a1*x+a2*x^2+����+apoly_n*x^poly_n
//n�����ݸ��� xy������ֵ poly_n�Ƕ���ʽ������
//����a0,a1,a2,����a[poly_n]���������У���matlab�෴
void polyfit(double a[], int n, double x[], double y[], int poly_n)
{
	int i, j;
	double *tempx, *tempy, *sumxx, *sumxy, *ata;
	void gauss_solve(int n, double A[], double x[], double b[]);
	tempx = calloc(n, sizeof(double));
	sumxx = calloc(poly_n * 2 + 1, sizeof(double));
	tempy = calloc(n, sizeof(double));
	sumxy = calloc(poly_n + 1, sizeof(double));
	ata = calloc((poly_n + 1)*(poly_n + 1), sizeof(double));
	for (i = 0; i<n; i++)
	{
		tempx[i] = 1;
		tempy[i] = y[i];
	}
	for (i = 0; i<2 * poly_n + 1; i++)
	for (sumxx[i] = 0, j = 0; j<n; j++)
	{
		sumxx[i] += tempx[j];
		tempx[j] *= x[j];
	}
	for (i = 0; i<poly_n + 1; i++)
	for (sumxy[i] = 0, j = 0; j<n; j++)
	{
		sumxy[i] += tempy[j];
		tempy[j] *= x[j];
	}
	for (i = 0; i<poly_n + 1; i++)
	for (j = 0; j<poly_n + 1; j++)
		ata[i*(poly_n + 1) + j] = sumxx[i + j];
	gauss_solve(poly_n + 1, ata, a, sumxy);
	free(tempx);
	free(sumxx);
	free(tempy);
	free(sumxy);
	free(ata);
}
void gauss_solve(int n, double A[], double x[], double b[])
{
	int i, j, k, r;
	double max;
	for (k = 0; k<n - 1; k++)
	{
		max = fabs(A[k*n + k]); /*find maxmum*/
		r = k;
		for (i = k + 1; i<n - 1; i++)
		if (max<fabs(A[i*n + i]))
		{
			max = fabs(A[i*n + i]);
			r = i;
		}
		if (r != k)
		for (i = 0; i<n; i++) /*change array:A[k]&A[r] */
		{
			max = A[k*n + i];
			A[k*n + i] = A[r*n + i];
			A[r*n + i] = max;
		}
		max = b[k]; /*change array:b[k]&b[r] */
		b[k] = b[r];
		b[r] = max;
		for (i = k + 1; i<n; i++)
		{
			for (j = k + 1; j<n; j++)
				A[i*n + j] -= A[i*n + k] * A[k*n + j] / A[k*n + k];
			b[i] -= A[i*n + k] * b[k] / A[k*n + k];
		}
	}
	for (i = n - 1; i >= 0; x[i] /= A[i*n + i], i--)
	for (j = i + 1, x[i] = b[i]; j<n; j++)
		x[i] -= A[i*n + j] * x[j];
}