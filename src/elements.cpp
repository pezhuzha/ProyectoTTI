#include "../include/elements.h"
#include "../include/SAT_Const.h"
#include <cmath>

    /**
     * @file elements.cpp
     * @brief El archivo contiene las implementaciones de elements.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	tuple<double,double,double,double,double,double,double> elements (Matrix y){ 
	 double pi2 = 2*pi,p=0.0, ap=0.0, e=0.0, i=0.0, Omega=0.0, omega=0.0,magh,eCosE,nu,E,M,e2,R,eSinE,H,u,a;
	 Matrix h;
	 
	 Matrix r = extract_vector(y,1,3);                                        
	 Matrix v = extract_vector(y,4,6);                                        
	 
	 h = cross(r,v);                                    
	 magh = norm(h);
	 p = magh*magh/GM_Earth;
	 H = norm(h);
	 
	 Omega = atan2 ( h(1), -h(2) );                     
	 Omega = fmod(Omega,pi2);
	 if(Omega<0){
	 	Omega+=pi2;
	 }
	 i     = atan2 ( sqrt(h(1)*h(1)+h(2)*h(2)), h(3) ); 
	 u     = atan2 ( r(3)*H, -r(1)*h(2)+r(2)*h(1) );    
	 
	 R  = norm(r);                                      
	 
	 a = 1/(2/R-dot(v,v)/GM_Earth);               
	 
	 eCosE = 1-R/a;                                     
	 eSinE = dot(r,v)/sqrt(GM_Earth*a);           
	 
	 e2 = eCosE*eCosE +eSinE*eSinE;
	 e  = sqrt(e2);                                     
	 E  = atan2(eSinE,eCosE);                           
	 
	 M  = fmod(E-eSinE,pi2);                             
	 
	 nu = atan2(sqrt(1.0-e2)*eSinE, eCosE-e2);          
	 
	 omega = fmod(u-nu,pi2);                             
	 
	 if(omega<0){
	 	omega+=pi2;
	 }
	 
		return tie(p, a, e, i, Omega, omega, M);
	}