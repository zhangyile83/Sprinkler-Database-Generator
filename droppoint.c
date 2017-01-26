#include <math.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "constant.h"
#define POINT_NUM 51
#define T_END 10
#define POLYN 4
int drop_without_wind(double t, const double y[], double f[], void *params);
int drop_with_wind(double t, const double y[], double f[], void *params);
double trapz(double t[], double x[], int length);
void polyfit(double a[], int n, double x[], double y[], int poly_n);
void roots(double a[], int poly_n, double z[]);
void droppoint(double *x, double *y, double v, double gama, double alpha, double beta, double diameter, double windspeed){
	
	double vx = v*cos((double)alpha*PI / 180)*cos((double)gama*PI / 180);
	double vy = v*cos((double)alpha*PI / 180)*sin((double)gama*PI / 180);
	double vz = v*sin((double)alpha*PI / 180);
	double t_i[POINT_NUM], vx_i[POINT_NUM], vy_i[POINT_NUM], vz_i[POINT_NUM];

	//ode start
	gsl_odeiv2_system sys = { drop_without_wind, NULL, 3, &diameter };
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);		
	double result[3];
	result[0] = vx; result[1] = vy; result[2] = vz;
	double tt_tmp=0;
	for (int i = 0; i < POINT_NUM; i++)
	{
		double ti = (double)i * T_END / (POINT_NUM - 1);
		int status = gsl_odeiv2_driver_apply(d, &tt_tmp, ti, result);
		if (status != GSL_SUCCESS)
			exit(0);//error		
		t_i[i] = tt_tmp;
		vz_i[i] = result[2];
	}
	gsl_odeiv2_driver_free(d);
	//ode end
	
	double z[POINT_NUM - 1];
	for (int i = 0; i < POINT_NUM - 1; i++)
		z[i] = trapz(t_i, vz_i, i+2);	
	double a[5];
	double roots_tmp[8];
	double drop_time=0;
	polyfit(a, POINT_NUM - 1, t_i+1, z, POLYN);
	roots(a, POLYN, roots_tmp);
	for (int i = 0; i < 4; i++){
		if (roots_tmp[2 * i + 1] == 0 && roots_tmp[2 * i]>0)
			drop_time = roots_tmp[2 * i];
	}
	if (drop_time == 0){
		exit(0);//error
	}
	//ode start	
	double para[3];
	para[0] = diameter;
	para[1] = windspeed;
	para[2] = beta;
	gsl_odeiv2_system sys2 = { drop_with_wind, NULL, 3, &para };
	gsl_odeiv2_driver * d2 = gsl_odeiv2_driver_alloc_y_new(&sys2, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);	
	result[0] = vx; result[1] = vy; result[2] = vz;
	tt_tmp = 0;
	for (int i = 0; i < POINT_NUM; i++)
	{
		double ti = (double)i * drop_time / (POINT_NUM - 1);
		int status = gsl_odeiv2_driver_apply(d2, &tt_tmp, ti, result);
		if (status != GSL_SUCCESS)
			exit(0);//error		
		t_i[i] = tt_tmp;
		vx_i[i] = result[0];
		vy_i[i] = result[1];
	}
	
	*x = trapz(t_i, vx_i, POINT_NUM);	//result[0];
	*y = trapz(t_i, vy_i, POINT_NUM);	//result[1];
	gsl_odeiv2_driver_free(d2);
	//ode end

	return;
}

int drop_without_wind(double t, const double y[], double f[], void *params)
{
	(void)(t); /* avoid unused parameter warning */
	double diameter = *(double *)params;
	double g = 9.8; // gravitation constant
	double fai = 0.18; // fai is a constant related to the Reynolds number
	double pa = 1.29; // pa is the atmospheric density
	double pw = 1000; // pw is the water density
	double m = 3.141592 / 6 * pw*diameter*diameter*diameter; // mass of the droplet????
	double k = fai*pa*diameter*diameter; // friction constant
	f[0] = -(k / m) * sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]) * y[0];
	f[1] = -(k / m) * sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]) * y[1];
	f[2] = -(k / m) * sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]) * y[2] - g;
	return GSL_SUCCESS;
}

int drop_with_wind(double t, const double y[], double f[], void *params)
{
	(void)(t); /* avoid unused parameter warning */
	double diameter = ((double *)params)[0];
	double windspeed = ((double *)params)[1];
	double beta = ((double *)params)[2];
	double g = 9.8; // gravitation constant
	double fai = 0.18; // fai is a constant related to the Reynolds number
	double pa = 1.29; // pa is the atmospheric density
	double pw = 1000; // pw is the water density
	double m = 3.141592 / 6 * pw*diameter*diameter*diameter; // mass of the droplet????
	double k = fai*pa*diameter*diameter; // friction constant
	double wx = windspeed * cos((double)beta*PI / 180);
	double wy = windspeed * sin((double)beta*PI / 180);
	f[0] = -(k / m) * sqrt((y[0] - wx) * (y[0] - wx) + (y[1] - wy) * (y[1] - wy) + y[2] * y[2]) * (y[0] - wx);	
	f[1] = -(k / m) * sqrt((y[0] - wx) * (y[0] - wx) + (y[1] - wy) * (y[1] - wy) + y[2] * y[2]) * (y[1] - wy);	
	f[2] = -(k / m) * sqrt((y[0] - wx) * (y[0] - wx) + (y[1] - wy) * (y[1] - wy) + y[2] * y[2]) * y[2] - g;	
	return GSL_SUCCESS;
}


double trapz(double t[], double x[], int length)
{
	double sum = 0;
	int i;
	for (i = 0; i<length - 1; i++)
	{
		sum += (t[i + 1] - t[i])*(x[i] + x[i + 1]) / 2;
	}
	return sum;//返回积分值
}