#include "../include/R_y.h"
#include "../include/matrix.h"

    /**
     * @file R_y.cpp
     * @brief El archivo contiene las implementaciones de R_y.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	Matrix& R_y(double angle){
		double C = cos(angle);
		double S = sin(angle);
		Matrix *rotmat=&zeros(3,3);

		(*rotmat)(1,1) =   C;  (*rotmat)(1,2) = 0.0;  (*rotmat)(1,3) = -1.0*S;
		(*rotmat)(2,1) = 0.0;  (*rotmat)(2,2) = 1.0;  (*rotmat)(2,3) =    0.0;
		(*rotmat)(3,1) =   S;  (*rotmat)(3,2) = 0.0;  (*rotmat)(3,3) =      C;
		return *rotmat;
	}