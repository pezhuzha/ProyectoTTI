#ifndef _VarEqn_
#define _VarEqn_
#include "matrix.h"

using namespace std;

    /**
     * @file VarEqn.h
     * @brief El archivo contiene la funcion VarEqn
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * Computes the variational equations, i.e. the derivative of the state vector and the state transition matrix
     * @param x           Time since epoch in [s]
     * @param yPhi        (6+36)-dim vector comprising the state vector (y) and the state transition matrix (Phi) in column wise storage order
     * @return yPhip       Derivative of yPhi
     */
	Matrix& VarEqn(double x, Matrix yPhi);
#endif


