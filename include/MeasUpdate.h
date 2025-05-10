#ifndef _MeasUpdate_
#define _MeasUpdate_
#include "matrix.h"

using namespace std;

    /**
     * @file MeasUpdate.h
     * @brief El archivo contiene la funcion MeasUpdate
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * @param x Matrix n*1
     * @param z double
     * @param g double
     * @param s double
     * @param G Matrix 1*n
     * @param P Matrix n*n
     * @param n int
     * @return tupla de 3 Matrix
     */
	tuple<Matrix&,Matrix&,Matrix&> MeasUpdate(Matrix x, double z,double g,double s,Matrix G,Matrix P, int n);
#endif


