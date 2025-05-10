#ifndef _Accel_
#define _Accel_
#include "matrix.h"

using namespace std;

    /**
     * @file Accel.h
     * @brief El archivo contiene la funcion Accel
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * Computes the acceleration of an Earth orbiting satellite due to 
     *    - the Earth's harmonic gravity field, 
     *    - the gravitational perturbations of the Sun and Moon
     *    - the solar radiation pressure and
     *    - the atmospheric drag
     * @param Mjd_TT      Terrestrial Time (Modified Julian Date)
     * @param Y           Satellite state vector in the ICRF/EME2000 system
     * @return dY           Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
     */
	Matrix& Accel(double x,Matrix Y);
#endif


