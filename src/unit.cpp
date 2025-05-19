#include "../include/unit.h"

    /**
     * @file unit.cpp
     * @brief El archivo contiene las implementaciones de unit.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	Matrix& unit(Matrix vec){
		int i=0;
		double small = 0.000001;
		double magv = norm(vec);
		Matrix &outvec=zeros(3);
		if ( magv > small ){
		    for (i=1;i<=3;i++){
		        outvec(i)= vec(i)/magv;
		    }
		}
		else{
		    for (i=1;i<=3;i++){
		    	outvec(i)= 0.0;
		    }
		}

		return outvec;
	}