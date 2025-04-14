#include "../include/sign_.h"
#include <cmath>

    /**
     * @file sign_.cpp
     * @brief El archivo contiene las implementaciones de sign_.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	double sign_(double a, double b){
	if (b>=0.0){
	    return abs(a);}
	else{
	    return -abs(a);}
	}