#include"constant.h"
double calc_error(double target_x, double target_y, double p1, double p2, double v0, double gama, double alpha, double beta, double diameter, double wind_speed);
void anti_wind_parameter(double *output_p1, double *output_p2, double v0, double target_x, double target_y, double alpha, double beta, double gama, double diameter, double wind_speed){
	double scale_p1 = 0.05;
	double scale_p2 = 0.002;
	static double p1 = 1;
	static double p2 = 1;
	static int adjust_direction_p1 = INCREASE;
	static int adjust_direction_p2 = INCREASE;
	double error_min = calc_error(target_x, target_y, p1, p2, v0, gama, alpha, beta, diameter, wind_speed);	

	do{
		//adjust p1
		int better = 0;
		double error_tmp;
		if ((error_tmp = calc_error(target_x, target_y, p1 + adjust_direction_p1*scale_p1, p2, v0, gama, alpha, beta, diameter, wind_speed)) < error_min)//%如果这个方向，则一直照这个方向运行直到极点
			better = 1;
		else if ((error_tmp = calc_error(target_x, target_y, p1 - adjust_direction_p1*scale_p1, p2, v0, gama, alpha, beta, diameter, wind_speed)) < error_min){
			adjust_direction_p1 = -adjust_direction_p1;
			better = 1;
		}
		if (better == 1){
			do {
				p1 = p1 + adjust_direction_p1*scale_p1;
				error_min = error_tmp;
				error_tmp = calc_error(target_x, target_y, p1 + adjust_direction_p1*scale_p1, p2, v0, gama, alpha, beta, diameter, wind_speed);
			} while (error_tmp<error_min);
		}
		//adjust p2
		better = 0;		
		if ((error_tmp = calc_error(target_x, target_y, p1, p2 + adjust_direction_p2*scale_p2, v0, gama, alpha, beta, diameter, wind_speed)) < error_min)//%如果这个方向，则一直照这个方向运行直到极点
			better = 1;
		else if ((error_tmp = calc_error(target_x, target_y, p1, p2 - adjust_direction_p2*scale_p2, v0, gama, alpha, beta, diameter, wind_speed)) < error_min){
			adjust_direction_p2 = -adjust_direction_p2;
			better = 1;
		}
		if (better == 1){
			do {
				p2 = p2 + adjust_direction_p2*scale_p2;
				error_min = error_tmp;
				error_tmp = calc_error(target_x, target_y, p1, p2 + adjust_direction_p2*scale_p2, v0, gama, alpha, beta, diameter, wind_speed);
			} while (error_tmp<error_min);
		}
		scale_p1 = scale_p1 / 2;
		scale_p2 = scale_p2 / 2;
	} while (error_min > THRESHOLD);
	*output_p1 = p1;
	*output_p2 = p2;
}