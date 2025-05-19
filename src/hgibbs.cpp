#include "../include/hgibbs.h"
#include "../include/unit.h"
#include "../include/angl.h"
#include "../include/SAT_Const.h"

    /**
     * @file hgibbs.cpp
     * @brief El archivo contiene las implementaciones de hgibbs.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	tuple <Matrix&,double,double,double,string>hgibbs(Matrix r1,Matrix r2,Matrix r3,double Mjd1, double Mjd2,double Mjd3){
		string error =  "ok";
		double theta = 0.0;
		double theta1= 0.0;
		double magr1 = norm( r1 );
		double magr2 = norm( r2 );
		double magr3 = norm( r3 );
		Matrix &v2=zeros(3);
		int i=0;

		double tolangle= 0.01745329251994;
		double dt21= (Mjd2-Mjd1)*86400.0;
		double dt31= (Mjd3-Mjd1)*86400.0;
		double dt32= (Mjd3-Mjd2)*86400.0;

		Matrix p = cross( r2,r3 );
		Matrix pn = unit( p );
		Matrix r1n = unit( r1 );
		double copa=  asin( dot( pn,r1n ) );

		if ( abs( dot(r1n,pn) ) > 0.017452406 ){
		    error= "not coplanar";
		}
		

		theta  = angl( r1,r2 );
		theta1 = angl( r2,r3 );

		if ( (theta > tolangle) | (theta1 > tolangle) )  {
			error= "angl > 1ø";
		}

		double term1= -dt32*( 1.0/(dt21*dt31) + GM_Earth/(12.0*magr1*magr1*magr1) );
		double term2= (dt32-dt21)*( 1.0/(dt21*dt32) + GM_Earth/(12.0*magr2*magr2*magr2) );
		double term3=  dt21*( 1.0/(dt32*dt31) + GM_Earth/(12.0*magr3*magr3*magr3) );

		v2 =  r1*term1 + r2*term2 + r3*term3 ;

		return tie(v2, theta,theta1,copa, error);
	}