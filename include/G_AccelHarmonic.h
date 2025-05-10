#ifndef _G_AccelHarmonic_
#define _G_AccelHarmonic_
#include "matrix.h"

using namespace std;

    /**
     * @file G_AccelHarmonic.h
     * @brief El archivo contiene la funcion G_AccelHarmonic
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * @param r           Satellite position vector in the true-of-date system
     * @param U           Transformation matrix to body-fixed syste
     * @param n           Gravity model degree
     * @param m             Gravity model order
     * @return G            Gradient (G=da/dr) in the true-of-date system
     */
	Matrix& G_AccelHarmonic( Matrix r,Matrix U, int n_max, int m_max );
#endif


