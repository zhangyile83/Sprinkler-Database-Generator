double angle_adjuster(double old_angle){
	while (old_angle >= 360)
		old_angle -= 360;
	while (old_angle < 0)
		old_angle += 360;
	return old_angle;
}