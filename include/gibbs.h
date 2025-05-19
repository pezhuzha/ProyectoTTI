#ifndef _gibbs_
#define _gibbs_
#include "matrix.h"
#include <tuple>
#include <cstring>
using namespace std;

    /**
     * @file gibbs.h
     * @brief El archivo contiene la funcion gibbs
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * @param r1          - ijk position vector #1         m
     * @param r2          - ijk position vector #2         m
     * @param r3          - ijk position vector #3         m
     * @return tuple[
     * v2          - ijk velocity vector for r2     m/s,
     * theta       - angl between vectors           rad,
     * theta1       double
     * copa         double
     * error       - flag indicating success        'ok',...
     * ]
     */
	tuple <Matrix&,double,double,double,string>gibbs(Matrix r1,Matrix r2,Matrix r3);
#endif


