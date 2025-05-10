#ifndef _GHAMatrix_
#define _GHAMatrix_
#include "matrix.h"

using namespace std;

    /**
     * @file GHAMatrix.h
     * @brief El archivo contiene la funcion GHAMatrix
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * Transformation from true equator and equinox to Earth equator and Greenwich meridian system 
     * @param Mjd_UT1   Modified Julian Date UT1
     * @return GHAmat    Greenwich Hour Angle matrix
     */
	Matrix& GHAMatrix (double Mjd_UT1);
#endif


