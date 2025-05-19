#ifndef _elements_
#define _elements_
#include "matrix.h"
#include <tuple>

using namespace std;

    /**
     * @file elements.h
     * @brief El archivo contiene la funcion elements
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * Computes the osculating Keplerian elements from the satellite state vector for elliptic orbits
     * The function cannot be used with state vectors describing a circular or non-inclined orbit.
     *  @param y        State vector (x,y,z,vx,vy,vz)
     * @return tuple<
     *  p        semilatus rectum [m],
     *  a        Semimajor axis ,
     *  e        Eccentricity ,
     *  i        Inclination [rad],
     *  Omega    Longitude of the ascending node [rad],
     *  omega    Argument of pericenter [rad],
     *  M        Mean anomaly [rad]
     *  >
     * @return Matrix resultado
     */
	tuple<double,double,double,double,double,double,double> elements (Matrix y);
#endif
