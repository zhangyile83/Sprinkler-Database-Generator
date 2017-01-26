#include <math.h>
void droppoint(double *x, double *y, double v, double gama, double alpha, double beta, double diameter, double windspeed);
double calc_error(double target_x, double target_y, double p1, double p2, double v0, double gama, double alpha, double beta, double diameter, double wind_speed){
	double drop_x, drop_y;
	droppoint(&drop_x, &drop_y, v0*p1, gama*p2, alpha, beta, diameter, wind_speed);
	double error = pow((target_x - drop_x) *(target_x - drop_x) + (target_y - drop_y) *(target_y - drop_y), 0.5);
	return error;
}