#ifndef _gast_
#define _gast_

using namespace std;

    /**
     * @file gast.h
     * @brief El archivo contiene la funcion gast
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * Greenwich Apparent Sidereal Time
     * @param Mjd_UT1   Modified Julian Date UT1
     * @return gstime    GAST in [rad]
     */
	double gast(double Mjd_UT1);
#endif


