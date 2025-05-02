#include "../include/LTC.h"
#include "../include/R_y.h"
#include "../include/R_z.h"
    /**
     * @file LTC.cpp
     * @brief El archivo contiene las implementaciones de LTC.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	Matrix& LTC(double lon,double lat){
		Matrix &M= R_y(-1.0*lat)*R_z(lon);
		double Aux;
	for (int j=1;j<=3;j++){
	    Aux=M(1,j);
	    M(1,j)=M(2,j);
	    M(2,j)=M(3,j);
	    M(3,j)= Aux;}
    return M;

	}