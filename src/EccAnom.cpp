#include "../include/EccAnom.h"
#include <cmath>
#include "../include/SAT_Const.h"
#include <iostream>
using namespace std;
    /**
     * @file EccAnom.cpp
     * @brief El archivo contiene las implementaciones de EccAnom.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	double EccAnom (double M,double e){

	double maxit = 15;
	int i = 1;

	M = fmod(M, 2.0*pi);
	double E;

	if (e<0.8){
		E = M; }
	else{
		E = pi;}

	double f = E - e*sin(E) - M;
	E = E - f / ( 1.0 - e*cos(E) );

	while (abs(f) > 1e2*eps)   {
		f = E - e*sin(E) - M;
		E = E - f / ( 1.0 - e*cos(E) );
		i = i+1;
		if (i==maxit){
			cerr<< "convergence problems in EccAnom";
		}
	}
	return E;
	}