#ifndef _AzElPa_
#define _AzElPa_
#include <tuple>
#include "matrix.h"
using namespace std;

    /**
     * @file AzElPa.h
     * @brief El archivo contiene la funcion AzElPa
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * Computes azimuth, elevation and partials from local tangent coordinates s  
     * @param s Matrix 1x3 Topocentric local tangent coordinates (East-North-Zenith frame)
     * @return tuple<A,E,dAds,dEds>      Azimuth [rad],Elevation [rad],Partials of azimuth w.r.t. s,Partials of elevation w.r.t. s
     */
	tuple<double,double,Matrix,Matrix> AzElPa(Matrix s);
#endif


