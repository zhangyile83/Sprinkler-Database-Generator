#define PRECISION_NUM 600
#define MAX_V 150
#include <math.h>
void droppoint(double *x, double *y, double v, double gama, double alpha, double beta, double diameter, double windspeed);
double characteristic(double alpha, double target_distance, double d){	
	static double distance[PRECISION_NUM];
	static double initialized_flag = 0;
	if (initialized_flag == 0){
		initialized_flag = 1;
		double x, y;
		for (int i = 0; i < PRECISION_NUM; i++){
			droppoint(&x, &y, (double)i*MAX_V / PRECISION_NUM, 360.0, alpha, 0.0, d, 0.0);
			distance[i] = x;
		}
	}
	
	double last = fabs(distance[0] - target_distance);
	for (int i = 1; i < PRECISION_NUM; i++){
		if (last>fabs(distance[i] - target_distance))
			last = fabs(distance[i] - target_distance);
		else
			return (double)(i-1)*MAX_V / PRECISION_NUM;
	}
	return 0;
}