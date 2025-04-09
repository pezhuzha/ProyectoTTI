#ifndef _Cheb3D_
#define _Cheb3D_
#include "matrix.h"

using namespace std;

    /**
     * @file Cheb3D.h
     * @brief El archivo contiene la funcion Cheb3D
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
	* @param     t       time
	* @param     N       Number of coefficients
	* @param     Ta      Begin interval
	* @param     Tb      End interval
	* @param     Cx      Coefficients of Chebyshev polyomial (x-coordinate)
	* @param     Cy      Coefficients of Chebyshev polyomial (y-coordinate)
	* @param     Cz      Coefficients of Chebyshev polyomial (z-coordinate)
     * @return Chebyshev approximation of 3-dimensional vectors
     */
	Matrix&	Cheb3D( double t, int N, double Ta, double Tb, Matrix& Cx,Matrix& Cy,Matrix& Cz);
#endif


