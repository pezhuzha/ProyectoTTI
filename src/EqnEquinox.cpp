#include "../include/EqnEquinox.h"
#include "../include/NutAngles.h"
#include "../include/MeanObliquity.h"
#include <cmath>


    /**
     * @file EqnEquinox.cpp
     * @brief El archivo contiene las implementaciones de EqnEquinox.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	double EqnEquinox (double Mjd_TT){

	auto [dpsi, deps] = NutAngles (Mjd_TT);

	return dpsi * cos ( MeanObliquity(Mjd_TT) );
}

