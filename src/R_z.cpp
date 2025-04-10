#include "../include/R_z.h"
#include "../include/matrix.h"

    /**
     * @file R_z.cpp
     * @brief El archivo contiene las implementaciones de R_z.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	Matrix& R_z(double angle){
		double C = cos(angle);
		double S = sin(angle);
		Matrix *rotmat=&zeros(3,3);
		
		(*rotmat)(1,1) =      C;  (*rotmat)(1,2) =   S;  (*rotmat)(1,3) = 0.0;
		(*rotmat)(2,1) = -1.0*S;  (*rotmat)(2,2) =   C;  (*rotmat)(2,3) = 0.0;
		(*rotmat)(3,1) =    0.0;  (*rotmat)(3,2) = 0.0;  (*rotmat)(3,3) = 1.0;
		return *rotmat;
	}