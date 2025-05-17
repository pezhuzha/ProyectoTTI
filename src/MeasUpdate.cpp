#include "../include/MeasUpdate.h"
    /**
     * @file MeasUpdate.cpp
     * @brief El archivo contiene las implementaciones de MeasUpdate.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	tuple<Matrix&,Matrix&,Matrix&> MeasUpdate(Matrix x, double z,double g,double s,Matrix G,Matrix P, int n){

			if(x.n_row<x.n_column){
				x=transpose(x);
			}
		Matrix Inv_W(1);Inv_W(1) = s*s;

		Matrix &K = P*transpose(G)*inv(Inv_W+G*P*transpose(G));

		Matrix &nx = x + K*(z-g);

		Matrix &nP = (eye(n)-K*G)*P;

		return tie(K, nx, nP);

	}