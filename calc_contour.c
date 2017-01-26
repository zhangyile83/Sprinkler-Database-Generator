#include"constant.h"
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
double calc_intersection(double angle, double x1, double y1, double x0, double y0);
int on_a_line_test(double x1, double y1, double x2, double y2);
void calc_contour(double contour_distance[], double *angle_start, double *angle_end, double vertex_x[], double vertex_y[], int vertex_number, double angle_step_length){
	double *angle;
	double tmpi;
	double tmpd;
	angle = (double*)malloc(sizeof(double)*vertex_number);
	for (int i = 0; i < vertex_number; i++){		
		angle[i] = (atan2(vertex_y[i], vertex_x[i]) * 180 / PI);
		if (angle[i] < 0)
			angle[i] += 360;
	}
	for (int i = 0; i < vertex_number-1;i++)
	for (int j = 0; j < vertex_number-1-i; j++){
		if (angle[j]>angle[j + 1]){
			tmpi = angle[j];
			angle[j] = angle[j + 1];
			angle[j + 1] = tmpi;
			tmpd = vertex_x[j];
			vertex_x[j] = vertex_x[j + 1];
			vertex_x[j + 1] = tmpd;
			tmpd = vertex_y[j];
			vertex_y[j] = vertex_y[j + 1];
			vertex_y[j + 1] = tmpd;
		}			
	}
	//judge if it's on a vertex
	int on_a_vertex = 0;
	int on_the_line = 0;
	for (int i = 0; i < vertex_number; i++)
	if (vertex_x[i] == 0 && vertex_y[i] == 0){
		on_a_vertex = 1;
		on_the_line = 1;
		vertex_number--;
		for (int j = i; j < vertex_number; j++){
			vertex_x[j] = vertex_x[j + 1];
			vertex_y[j] = vertex_y[j + 1];
			angle[j] = angle[j + 1];
		}
		break;
	}
	if (on_a_vertex == 0){
		for (int i = 0; i < vertex_number - 1; i++)
		if (on_a_line_test(vertex_x[i], vertex_y[i], vertex_x[i + 1], vertex_y[i + 1])){
			on_the_line = 1;
			break;
		}
		if (on_a_line_test(vertex_x[vertex_number - 1], vertex_y[vertex_number - 1], vertex_x[0], vertex_y[0])){
			on_the_line = 1;
		}
	}
	if (on_the_line == 0){
		*angle_start = 0;
		*angle_end = 359;
	}
	else{
		*angle_start = angle[0];
		*angle_end = angle[vertex_number - 1];
	}
	int count_tmp = 0;
	double angle_tmp = 0;
	if (on_the_line == 0){
		for (; angle_tmp < angle[0]; angle_tmp += angle_step_length){
			contour_distance[count_tmp] = calc_intersection(angle_tmp, vertex_x[vertex_number - 1], vertex_y[vertex_number - 1], vertex_x[0], vertex_y[0]);
			count_tmp++;
		}
		for (int vertex_index = 0; vertex_index < vertex_number - 1; vertex_index++){
			for (; angle_tmp < angle[vertex_index + 1]; angle_tmp += angle_step_length){
				contour_distance[count_tmp] = calc_intersection(angle_tmp, vertex_x[vertex_index], vertex_y[vertex_index], vertex_x[vertex_index + 1], vertex_y[vertex_index + 1]);
				count_tmp++;
			}
		}
		for (; angle_tmp <=359; angle_tmp += angle_step_length){
			contour_distance[count_tmp] = calc_intersection(angle_tmp, vertex_x[vertex_number - 1], vertex_y[vertex_number - 1], vertex_x[0], vertex_y[0]);
			count_tmp++;
		}
	}
	else{
		for (int vertex_index = 0; vertex_index < vertex_number - 1; vertex_index++){
			for (angle_tmp = angle[vertex_index]; angle_tmp < angle[vertex_index + 1]; angle_tmp += angle_step_length){
				contour_distance[count_tmp] = calc_intersection(angle_tmp, vertex_x[vertex_index], vertex_y[vertex_index], vertex_x[vertex_index + 1], vertex_y[vertex_index + 1]);
				count_tmp++;
			}
		}
	}
	*angle_start += 360;
	*angle_end += 360;
}

int on_a_line_test(double x1, double y1, double x2, double y2){
	if (x1 != 0 && y1 != 0 && x2 != 0 && y2 != 0 && y1*x2 == y2*x1)
		return 1;
	if (x1 == 0 && x2 == 0)
		return 1;
	if (y1 == 0 && y2 == 0)
		return 1;
	return 0;
}

double calc_intersection(double angle, double x1, double y1, double x0, double y0){
	double k1, k2, x, y;
	if (x0 - x1 != 0){
		k1 = (y0 - y1) / (x0 - x1);
		k2 = tan((double)angle*PI / 180); //90¶È½Ç²»ÐÐ
		if (k2 - k1 == 0){
			exit(0);
		}
		x = (k1*x0 - k2 * 0 + 0 - y0) / (k1 - k2);
		y = y0 + (x - x0)*k1;
	}
	else{
		x = x0;
		y = tan((double)angle*PI / 180) * x;
	}
	return sqrt(x*x + y*y);
}