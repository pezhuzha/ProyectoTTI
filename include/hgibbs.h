#ifndef _hgibbs_
#define _hgibbs_
#include "matrix.h"
#include <tuple>
#include <cstring>
using namespace std;

    /**
     * @file hgibbs.h
     * @brief El archivo contiene la funcion hgibbs
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * @param r1          - ijk position vector #1         m
     * @param r2          - ijk position vector #2         m
     * @param r3          - ijk position vector #3         m
     * @param Mjd1        - julian date of 1st sighting    days from 4713 bc
     * @param Mjd2        - julian date of 2nd sighting    days from 4713 bc
     * @param Mjd3        - julian date of 3rd sighting    days from 4713 bc
     * @return tuple[
     * v2          - ijk velocity vector for r2     m/s,
     * theta       - angl between vectors           rad,
     * theta1       double
     * copa         double
     * error       - flag indicating success        'ok',...
     * ]
     */
	tuple <Matrix&,double,double,double,string>hgibbs(Matrix r1,Matrix r2,Matrix r3,double Mjd1, double Mjd2,double Mjd3);
#endif


