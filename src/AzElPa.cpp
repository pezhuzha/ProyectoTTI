#include "../include/AzElPa.h"
#include "../include/SAT_Const.h"
    /**
     * @file AzElPa.cpp
     * @brief El archivo contiene las implementaciones de AzElPa.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	tuple<double,double,Matrix&,Matrix&> AzElPa(Matrix s) {

		 double rho = sqrt(s(1)*s(1)+s(2)*s(2));

		// Angles
		double Az = atan2(s(1),s(2));

		if (Az<0.0) {
		    Az = Az+pi2;
		}

		double El = atan ( s(3) / rho );

		// Partials
		Matrix &dAds = zeros(3);
		dAds(1)=s(2)/(rho*rho);
		dAds(2)=-s(1)/(rho*rho);
		dAds(3)=0.0 ;
		Matrix &dEds= zeros(3);
		dEds(1)=-s(1)*s(3)/rho;
		dEds(2)=-s(2)*s(3)/rho;
		dEds(3)=rho;
		 dEds= dEds/ dot(s,s);
		return tie(Az, El, dAds, dEds);

	}