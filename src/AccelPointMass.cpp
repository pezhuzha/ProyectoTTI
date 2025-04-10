#include "../include/AccelPointMass.h"
#include "../include/matrix.h"
#include <cmath>

    /**
     * @file AccelPointMass.cpp
     * @brief El archivo contiene las implementaciones de AccelPointMass.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	Matrix&	AccelPointMass(Matrix& r,Matrix& s,double GM){
		Matrix d = r - s;
		Matrix *a=new Matrix(3);
		*a=( d/pow(norm(d),3) + s/(pow(norm(s),3)) * -GM);

		return *a;
	}