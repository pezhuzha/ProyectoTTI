#ifndef _IERS_
#define _IERS_
#include <tuple>
#include "matrix.h"
using namespace std;

    /**
     * @file IERS.h
     * @brief El archivo contiene la funcion IERS
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * IERS: Management of IERS time and polar motion data
     * @param eop matrix 4x13
     * @param Mjd_UTC double
     * @param interp char
     * @return tupla <x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC>
     */
	tuple<double,double,double,double,double,double,double,double,double> IERS(Matrix eop,double Mjd_UTC,char interp='n');
#endif


