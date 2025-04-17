#ifndef _NutAngles_
#define _NutAngles_
#include <tuple>
using namespace std;

    /**
     * @file NutAngles.h
     * @brief El archivo contiene la funcion NutAngles
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * Nutation in longitude and obliquity
     * @param Mjd_TT     Modified Julian Date (Terrestrial Time)
     * @return tupla < dpsi,deps>  Nutation Angles
     */
	tuple<double,double> NutAngles (double Mjd_TT);
#endif


