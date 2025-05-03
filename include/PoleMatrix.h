#ifndef _PoleMatrix_
#define _PoleMatrix_
#include "matrix.h"

using namespace std;

    /**
     * @file PoleMatrix.h
     * @brief El archivo contiene la funcion PoleMatrix
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * Transformation from pseudo Earth-fixed to Earth-fixed coordinates for a given date
     * @param Pole coordinte(xp,yp)
     * @return PoleMat   Pole matrix
     */
	Matrix& PoleMatrix (double xp,double yp);
#endif


