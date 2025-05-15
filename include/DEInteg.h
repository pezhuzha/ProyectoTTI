#ifndef _DEInteg_
#define _DEInteg_
#include "matrix.h"

using namespace std;

    /**
     * @file DEInteg.h
     * @brief El archivo contiene la funcion DEInteg
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
	 * Numerical integration methods for ordinaray differential equations
	 *  This module provides implemenation of the variable order variable stepsize multistep method of Shampine & Gordon.
     * @param f funcion pasas double y Matrix devuelve Matrix
	 * @param t double
	 * @param tout double
	 * @param relerr double
	 * @param abserr double
	 * @param n_eqn int
	 * @param y Matrix
     * @return Matrix resultado
     */
	 Matrix& DEInteg(Matrix& f(double t,Matrix y),double t, double tout,double relerr,double abserr,int n_eqn,Matrix &y);
#endif


