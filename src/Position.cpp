#include "../include/Position.h"
#include "../include/SAT_Const.h"

    /**
     * @file Position.cpp
     * @brief El archivo contiene las implementaciones de Position.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	Matrix& Position(double lon,double lat,double h){
		double R_equ = R_Earth;
		double f     = f_Earth;

		double e2     = f*(2.0-f);
		double CosLat = cos(lat);
		double SinLat = sin(lat);

		double N = R_equ / sqrt(1.0-e2*SinLat*SinLat);

		Matrix *r=new Matrix(3);

		(*r)(1) =  (         N+h)*CosLat*cos(lon);
		(*r)(2) =  (         N+h)*CosLat*sin(lon);
		(*r)(3) =  ((1.0-e2)*N+h)*SinLat;

		return *r;
	}