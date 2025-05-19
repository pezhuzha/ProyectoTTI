#include "../include/Geodetic.h"
#include "../include/SAT_Const.h"
#include <iostream>
    /**
     * @file Geodetic.cpp
     * @brief El archivo contiene las implementaciones de Geodetic.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	tuple<double,double,double> Geodetic(Matrix r){
	double	R_equ = R_Earth;
	double 	f     = f_Earth;

	double epsRequ = eps*R_equ;        
	double e2      = f*(2.0-f);        

	double X = r(1);                   
	double Y = r(2);
	double Z = r(3);
	double rho2 = X*X + Y*Y;           
	double lon=0.0,lat=0.0,h=0.0,dZ=0.0,ZdZ=0.0,Nh=0.0,SinPhi=0.0,N=0.0,dZ_new=0.0;

	if (norm(r)==0.0){
		cerr<<"invalid input in Geodetic constructor\n";
		 lon = 0.0;
		 lat = 0.0;
		  h   = -R_Earth;
	}
	dZ = e2*Z;

	while(true){
		    ZdZ    =  Z + dZ;
		    Nh     =  sqrt ( rho2 + ZdZ*ZdZ ); 
		    SinPhi =  ZdZ / Nh;                    
		    N      =  R_equ / sqrt(1.0-e2*SinPhi*SinPhi);
		    dZ_new =  N*e2*SinPhi;
		    if ( abs(dZ-dZ_new) < epsRequ ){
		    	break;
		    		    }
		    dZ = dZ_new;
		}

	lon = atan2 ( Y, X );
	lat = atan2 ( ZdZ, sqrt(rho2) );
	h   = Nh - N;

		return tie(lon, lat, h);
	}