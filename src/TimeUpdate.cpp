#include "../include/TimeUpdate.h"

    /**
     * @file TimeUpdate.cpp
     * @brief El archivo contiene las implementaciones de TimeUpdate.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	Matrix TimeUpdate(Matrix P,Matrix Phi,Matrix Qdt){

		return Phi*P*transpose(Phi)+ Qdt;

	}

	Matrix TimeUpdate(Matrix P,Matrix Phi){
		Matrix Qdt=zeros(P.n_column,P.n_row);
		return Phi*P*transpose(Phi)+ Qdt;

	}