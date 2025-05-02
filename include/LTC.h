#ifndef _LTC_
#define _LTC_
#include "matrix.h"

using namespace std;

    /**
     * @file LTC.h
     * @brief El archivo contiene la funcion LTC
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * Transformation from Greenwich meridian system to local tangent coordinates
     * @param lon      -Geodetic East longitude [rad]
     * @param lat      -Geodetic latitude [rad]
     * @return M        -Rotation matrix from the Earth equator and Greenwich meridian to the local tangent (East-North-Zenith) coordinate system
     */
    Matrix& LTC(double lon,double lat);
#endif


