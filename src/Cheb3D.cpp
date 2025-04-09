#include "../include/Cheb3D.h"
#include "../include/matrix.h"

    /**
     * @file Cheb3D.cpp
     * @brief El archivo contiene las implementaciones de Cheb3D.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	Matrix& Cheb3D(t, N, Ta, Tb, Cx, Cy, Cz){
		
		if ( (t<Ta) || (Tb<t) ){
			cerr<<'ERROR: Time out of range in Cheb3D::Value\n';
			exit(EXIT_FAILURE);
			}
			
		double tau = (2*t-Ta-Tb)/(Tb-Ta);  

		Matrix f1 = zeros(1,3);
		Matrix f2 = zeros(1,3);
		Matrix old_f1= zeros(1,3);
		Matrix aux= zeros(1,3);
		for (i=N;i>2;i--){
			old_f1 = f1;
			aux(1)=Cx(i);
			aux(2)=Cy(i);
			aux(3)=Cz(i);
			f1 = 2*tau*f1-f2+aux;
			f2 = old_f1;
		}
		Matrix *ChebApp;
		*ChebApp= zeros(1,3);
		*ChebApp = tau*f1-f2+[Cx(1),Cy(1),Cz(1)];
		return *ChebApp;
	}