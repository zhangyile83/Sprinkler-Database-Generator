#include <stdio.h>
#include "constant.h"
#include <math.h>
void calc_contour(double contour_distance[], double *angle_start, double *angle_end, double vertex_x[], double vertex_y[], int vertex_number, double angle_step_length);
double angle_adjuster(double old_angle);
void anti_wind_parameter(double *output_p1, double *output_p2, double v0, double target_x, double target_y, double alpha, double beta, double gama, double diameter, double wind_speed);
double characteristic(double alpha, double target_distance, double d);
void main()
{
	////debug
	FILE *fp;
	fopen_s(&fp, "C:\\Users\\Yile\\Desktop\\Results.txt", "wb+");

	//user input
	double alpha = 30;  // The angle of elevation of the nozzle
	double MeanD = 0.01;  // Droplet mean diameter(should be the value defined in the range)	
	int vertex_number = 4; //Number of vertex
	double vertex_x[4] = { 22.86, 22.86, -22.86, -22.86 }; //Location of vertex
	double vertex_y[4] = { 0, 18.28, 18.28, 0 };
	double distance_proportion_start = (double)1 / 23, distance_proportion_step_length = (double)1 / 23, //target distance
		distance_proportion_end = 1,
		wind_speed_start = 0, wind_speed_step_length = 0.447, wind_speed_end = 6.7050, //wind direction
		beta_start = 0, beta_step_length = 5, beta_end = 359; //wind diretion
	double angle_step_length = 1;
	//calculation
	double p1, p2;
	int count = 0;
	double contour_distance[360];
	double angle_start, angle_end;
	calc_contour(contour_distance, &angle_start, &angle_end, vertex_x, vertex_y, vertex_number, angle_step_length);
	for (double distance_proportion = distance_proportion_start; distance_proportion <= distance_proportion_end; distance_proportion += distance_proportion_step_length){
		for (double wind_speed = wind_speed_start; wind_speed <= wind_speed_end; wind_speed += wind_speed_step_length){
			for (double beta = beta_start; beta <= beta_end; beta += beta_step_length){
				int index_gama = 0;
				for (double gama = angle_start; gama <= angle_end; gama += angle_step_length, index_gama++){
 					double distance = contour_distance[index_gama] * distance_proportion;
					double v0 = characteristic(alpha, distance, MeanD);
					anti_wind_parameter(&p1, &p2, v0, distance*cos((double)gama*PI / 180), distance*sin((double)gama*PI / 180), alpha, beta, gama, MeanD, wind_speed);
					//output					
					count++;
					printf("Case %d:\ndistance proportion=%f\twind speed=%f\nwind direction=%f\ttarget direction=%f\n", count, distance_proportion, wind_speed, beta, angle_adjuster(gama));
					printf("Solution %d:\nvelocity=%f\tdirection=%f\t\n\n", count, v0*p1, angle_adjuster(gama*p2));
					//debug
					//printf("%f\n", (double)count / 4795776);
					//fprintf(fp,"%f\t%f\n", v0*p1, angle_adjuster(gama*p2));					
					fprintf(fp, "%f\t%f\t%f\t%f\t%f\t%f\n", distance_proportion, wind_speed, beta, angle_adjuster(gama), v0*p1, angle_adjuster(gama*p2));
					if (count > 10){
						system("pause");
						fclose(fp);
						exit(0);
					}
				}
			}
		}
	}
	system("pause");
	fclose(fp);
}
//debug
//fclose(fp);

