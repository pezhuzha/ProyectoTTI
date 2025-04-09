#ifndef _AccelPointMass_
#define _AccelPointMass_
#include "matrix.h"

using namespace std;

    /**
     * @file AccelPointMass.h
     * @brief El archivo contiene la funcion AccelPointMass
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * @param r           Satellite position vector 
     * @param s           Point mass position vector
     * @param GM          Gravitational coefficient of point mass
     * @return Acceleration (a=d^2r/dt^2)
     */
	double	AccelPointMass(Matrix& r,Matrix& s,double GM);
#endif


