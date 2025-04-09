#include "../include/AccelPointMass.h"
#include "../include/matrix.h"

    /**
     * @file AccelPointMass.cpp
     * @brief El archivo contiene las implementaciones de AccelPointMass.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	double	AccelPointMass(Matrix& r,Matrix& s,double GM){
		
		double d = r - s;
		double a = -GM * ( d/(norm(d)**3) + s/(norm(s)**3) );

		return a;
	}