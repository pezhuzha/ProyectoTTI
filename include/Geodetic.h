#ifndef _Geodetic_
#define _Geodetic_
#include "matrix.h"
#include <tuple>

using namespace std;

    /**
     * @file Geodetic.h
     * @brief El archivo contiene la funcion Geodetic
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * geodetic coordinates (Longitude [rad], latitude [rad], altitude [m]) from given position vector (r [m])
     * @param r vector posicion
     * @return tuple<
     * lon          longitud,
     * lat          latitud,
     * h            altitud
     * >
     */
	tuple<double,double,double> Geodetic(Matrix r);
#endif


