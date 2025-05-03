#ifndef _NutMatrix_
#define _NutMatrix_
#include "matrix.h"

using namespace std;

    /**
     * @file NutMatrix.h
     * @brief El archivo contiene la funcion NutMatrix
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * Transformation from mean to true equator and equinox
     * @param Mjd_TT    Modified Julian Date (Terrestrial Time)
     * @return NutMat    Nutation matrix
     */
	Matrix& NutMatrix (double Mjd_TT);
#endif


