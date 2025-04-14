#include "../include/MeanObliquity.h"
#include "../include/matrix.h"
#include "../include/SAT_Const.h"

    /**
     * @file MeanObliquity.cpp
     * @brief El archivo contiene las implementaciones de MeanObliquity.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	 double MeanObliquity (double Mjd_TT){


		double T = (Mjd_TT-MJD_J2000)/36525;

		return Rad *( 84381.448/3600-(46.8150+(0.00059-0.001813*T)*T)*T/3600 );
	}