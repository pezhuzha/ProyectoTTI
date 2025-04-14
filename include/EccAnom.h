#ifndef _EccAnom_
#define _EccAnom_

using namespace std;

    /**
     * @file EccAnom.h
     * @brief El archivo contiene la funcion EccAnom
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
	 *Computes the eccentric anomaly for elliptic orbits
     * @param  M         Mean anomaly in [rad]
     * @param  e         Eccentricity of the orbit [0,1]
     * @return double resultado
     */
	double EccAnom (double M,double e);
#endif


