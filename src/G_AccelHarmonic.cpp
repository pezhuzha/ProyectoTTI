#include "../include/G_AccelHarmonic.h"
#include "../include/AccelHarmonic.h"

    /**
     * @file G_AccelHarmonic.cpp
     * @brief El archivo contiene las implementaciones de G_AccelHarmonic.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	Matrix& G_AccelHarmonic( Matrix r,Matrix U, int n_max, int m_max ){
		Matrix da;
	double d = 1.0;   

	Matrix &G = zeros(3,3);
	Matrix dr = zeros(3,1);

	for (int i=1;i<=3;i++){
		    dr = zeros(3,1);
		    dr(i) = d/2;

		    da = AccelHarmonic ( r+dr,U, n_max, m_max ) - 
		    AccelHarmonic ( r-dr,U, n_max, m_max );
		    G=assign_column(G,da/d,i) ;
		       } 
	return G;
	}