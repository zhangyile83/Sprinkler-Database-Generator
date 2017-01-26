#include <gsl/gsl_poly.h>
//求解0=a0+a1*x+a2*x^2+……+apoly_n*x^poly_n
//a[0],a[1],a[2],……a[poly_n]是系数，升序排列，跟matlab相反
void roots(double a[], int poly_n, double z[]){	
	gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc(poly_n+1);
	gsl_poly_complex_solve(a, poly_n + 1, w, z);
	gsl_poly_complex_workspace_free(w);
}