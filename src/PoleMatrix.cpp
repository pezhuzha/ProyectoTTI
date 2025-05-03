#include "../include/PoleMatrix.h"
#include "../include/R_x.h"
#include "../include/R_y.h"

    /**
     * @file PoleMatrix.cpp
     * @brief El archivo contiene las implementaciones de PoleMatrix.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	Matrix& PoleMatrix (double xp,double yp){
return R_y(-xp) * R_x(-yp);
}