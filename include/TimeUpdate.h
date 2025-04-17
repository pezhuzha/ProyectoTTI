#ifndef _TimeUpdate_
#define _TimeUpdate_
using namespace std;
#include "../include/matrix.h"

    /**
     * @file TimeUpdate.h
     * @brief El archivo contiene la funcion TimeUpdate
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * @param P Matrix
     * @param Phi Matrix
     * @param Qdt Matrix
     * @return Matrix
     */
	Matrix TimeUpdate(Matrix P,Matrix Phi,Matrix Qdt);

    /**
     * @param P Matrix
     * @param Phi Matrix
     * @return Matrix
     */
    Matrix TimeUpdate(Matrix P,Matrix Phi);
#endif


