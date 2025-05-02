#ifndef _JPL_Eph_DE430_
#define _JPL_Eph_DE430_
using namespace std;
#include"../include/matrix.h"

    /**
     * @file JPL_Eph_DE430.h
     * @brief El archivo contiene la funcion JPL_Eph_DE430
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * @param P Matrix
     * @param Phi Matrix
     * @param Qdt Matrix
     * @return Matrix
     */
	tuple<Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&> JPL_Eph_DE430(double Mjd_TDB);

#endif


