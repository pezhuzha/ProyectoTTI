#include "../include/Accel.h"
#include "../include/PrecMatrix.h"
#include "../include/NutMatrix.h"
#include "../include/IERS.h"
#include "../include/timediff.h"
#include "../include/PoleMatrix.h"
#include "../include/AccelHarmonic.h"
#include "../include/GHAMatrix.h"
#include "../include/JPL_Eph_DE430.h"
#include "../include/GLOBAL.h"
#include "../include/AccelPointMass.h"
#include "../include/Mjday_TDB.h"
#include "../include/SAT_Const.h"
    /**
     * @file Accel.cpp
     * @brief El archivo contiene las implementaciones de Accel.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
		Matrix& Accel(double x,Matrix Y){
			AuxParamLoad();
	auto [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,AuxParam.Mjd_UTC + x/86400,'l');
	auto [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
	double Mjd_UT1 = AuxParam.Mjd_UTC + x/86400 + UT1_UTC/86400;
	double Mjd_TT = AuxParam.Mjd_UTC + x/86400 + TT_UTC/86400;

	Matrix P = PrecMatrix(MJD_J2000,Mjd_TT);
	Matrix N = NutMatrix(Mjd_TT);
	Matrix T = N * P;
	Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

	double MJD_TDB = Mjday_TDB(Mjd_TT);
	auto [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,r_Sun] = JPL_Eph_DE430(MJD_TDB);

	// Acceleration due to harmonic gravity field
	Matrix a = AccelHarmonic(extract_vector(Y,1,3), E, AuxParam.n, AuxParam.m);

	// Luni-solar perturbations
	if (AuxParam.sun){
		    a = a + AccelPointMass(extract_vector(Y,1,3),r_Sun,GM_Sun);}

	if (AuxParam.moon){
		    a = a + AccelPointMass(extract_vector(Y,1,3),r_Moon,GM_Moon);}

	// Planetary perturbations
	if (AuxParam.planets){
		    a = a + AccelPointMass(extract_vector(Y,1,3),r_Mercury,GM_Mercury);
		    a = a + AccelPointMass(extract_vector(Y,1,3),r_Venus,GM_Venus);
		    a = a + AccelPointMass(extract_vector(Y,1,3),r_Mars,GM_Mars);
		    a = a + AccelPointMass(extract_vector(Y,1,3),r_Jupiter,GM_Jupiter);
		    a = a + AccelPointMass(extract_vector(Y,1,3),r_Saturn,GM_Saturn);
		    a = a + AccelPointMass(extract_vector(Y,1,3),r_Uranus,GM_Uranus);    
		    a = a + AccelPointMass(extract_vector(Y,1,3),r_Neptune,GM_Neptune);
		    a = a + AccelPointMass(extract_vector(Y,1,3),r_Pluto,GM_Pluto);}

	return union_vector(extract_vector(Y,1,3),a);

	}