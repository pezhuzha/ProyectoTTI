#ifndef _PrecMatrix_
#define _PrecMatrix_
#include "matrix.h"

using namespace std;

    /**
     * @file PrecMatrix.h
     * @brief El archivo contiene la funcion PrecMatrix
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * Precession transformation of equatorial coordinates
     * @param Mjd_1     Epoch given (Modified Julian Date TT)
     * @param MjD_2     Epoch to precess to (Modified Julian Date TT)
     * @return PrecMat   Precession transformation matrix
     */
	Matrix& PrecMatrix (double Mjd_1, double Mjd_2);
#endif


