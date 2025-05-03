#ifndef _gmst_
#define _gmst_

using namespace std;

    /**
     * @file gmst.h
     * @brief El archivo contiene la funcion gmst
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * Greenwich Mean Sidereal Time
     * @param Mjd_UT1    Modified Julian Date UT1
     * @return gmstime     GMST in [rad]
     */
	double gmst(double Mjd_UT1);
#endif


