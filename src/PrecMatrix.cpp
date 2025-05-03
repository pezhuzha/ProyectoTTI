#include "../include/PrecMatrix.h"
#include "../include/R_x.h"
#include "../include/R_y.h"
#include "../include/R_z.h"
#include "../include/SAT_Const.h"

    /**
     * @file PrecMatrix.cpp
     * @brief El archivo contiene las implementaciones de PrecMatrix.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	Matrix& PrecMatrix (double Mjd_1, double Mjd_2){
        double zeta,z,theta,dT,T;
        T  = (Mjd_1-MJD_J2000)/36525;
        dT = (Mjd_2-Mjd_1)/36525;

        
        zeta  =  ( (2306.2181+(1.39656-0.000139*T)*T)+ ((0.30188-0.000344*T)+0.017998*dT)*dT )*dT/Arcs;
        z     =  zeta + ( (0.79280+0.000411*T)+0.000205*dT)*dT*dT/Arcs;
        theta =  ( (2004.3109-(0.85330+0.000217*T)*T)-((0.42665+0.000217*T)+0.041833*dT)*dT )*dT/Arcs;

        return R_z(-z) * R_y(theta) * R_z(-zeta);
}