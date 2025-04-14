#ifndef _Position_
#define _Position_
#include "matrix.h"

using namespace std;

    /**
     * @file Position.h
     * @brief El archivo contiene la funcion Position
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * @param lon   longitude [rad]
     * @param lat   latitude [rad]
     * @param h     altitude [m]
     * @return Position vector (r [m]) from geodetic coordinates (Longitude [rad],latitude [rad], altitude [m])
     */
	Matrix& Position(double lon,double lat,double h);
#endif


