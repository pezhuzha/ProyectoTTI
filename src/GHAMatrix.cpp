#include "../include/GHAMatrix.h"
#include "../include/gast.h"
#include "../include/R_z.h"

    /**
     * @file GHAMatrix.cpp
     * @brief El archivo contiene las implementaciones de GHAMatrix.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	Matrix& GHAMatrix (double Mjd_UT1){
		return  R_z( gast(Mjd_UT1) );
	}