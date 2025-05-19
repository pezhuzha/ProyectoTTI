#include "../include/angl.h"
#include "../include/sign_.h"
#include <cmath>
    /**
     * @file angl.cpp
     * @brief El archivo contiene las implementaciones de angl.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	double angl ( Matrix vec1,Matrix vec2 ){
		
		double small     = 0.00000001;
		double undefined = 999999.1;
		double theta=0.0;
		double temp=0.0;
		double magv1 = norm(vec1);
		double magv2 = norm(vec2);

		if (magv1*magv2 > small*small){
		    temp= dot(vec1,vec2) / (magv1*magv2);
		    if (abs( temp ) > 1.0){
		        temp= sign_(temp,temp) * 1.0;
		    }
		    theta= acos( temp );}
		else{
			theta= undefined;
		}
		return theta;
	}