#ifndef _timediff_
#define _timediff_
#include <tuple>
using namespace std;

    /**
     * @file timediff.h
     * @brief El archivo contiene la funcion timediff
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * Time differences [s]
     * @param UT1_UTC double
     * @param TAI_UTC double
     * @return tupla <UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC>
     */
	tuple<double,double,double,double,double> timediff(double UT1_UTC,double TAI_UTC);
#endif


