#ifndef _AccelHarmonic_
#define _AccelHarmonic_
using namespace std;
#include "../include/matrix.h"

    /**
     * @file AccelHarmonic.h
     * @brief El archivo contiene la funcion AccelHarmonic
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	 
    /**
     * @param r           Satellite position vector in the inertial system  
     * @param E           Transformation matrix to body-fixed system  
     * @param n_max       Maximum degree  
     * @param m_max       Maximum order (m_max<=n_max; m_max=0 for zonals, only)
     * @return a           Acceleration (a=d^2r/dt^2)
     */
	double AccelHarmonic(Matrix r,Matrix E,int n_max,int m_max);
#endif


