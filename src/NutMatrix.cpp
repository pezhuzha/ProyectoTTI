#include "../include/NutMatrix.h"
#include "../include/MeanObliquity.h"
#include "../include/NutAngles.h"
#include "../include/R_x.h"
#include "../include/R_z.h"

    /**
     * @file NutMatrix.cpp
     * @brief El archivo contiene las implementaciones de NutMatrix.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
	Matrix& NutMatrix (double Mjd_TT){

double eps = MeanObliquity (Mjd_TT);

auto [dpsi, deps] = NutAngles (Mjd_TT);

return R_x(-eps-deps)*R_z(-dpsi)*R_x(+eps);
}