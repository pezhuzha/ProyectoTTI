#ifndef _Legendre_
#define _Legendre_
#include <tuple>
#include "matrix.h"
using namespace std;

    /**
     * @file Legendre.h
     * @brief El archivo contiene la funcion Legendre
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * Time differences [s]
     * @param n int
     * @param m int
     * @param fi double [rad]
     * @return tupla <pnm,dpnm>
     */
	tuple<Matrix&,Matrix&> Legendre(int n,int m,double fi);
#endif


