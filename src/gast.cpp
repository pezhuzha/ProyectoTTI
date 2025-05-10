#include "../include/gast.h"
#include "../include/gmst.h"
#include "../include/EqnEquinox.h"
#include "../include/SAT_Const.h"
#include <cmath>

    /**
     * @file gast.cpp
     * @brief El archivo contiene las implementaciones de gast.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	double gast(double Mjd_UT1){
		return fmod(gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), 2*pi );
	}