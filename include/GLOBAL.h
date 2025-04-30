    /**
     * @file GLOBAL.h
     * @brief El archivo contiene la funcion GLOBAL
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
#ifndef _GLOBAL_
#define _GLOBAL_
using namespace std;
#include "matrix.h"
#include <cmath>

extern Matrix eopdata;

    /**
     * Lee el archivo eop19620101.txt y recoge cada fila y lo asigna a eopdata
     * @param c    número de filas a recoger
     */
void eop19620101(int c);

extern Matrix Cnm;
extern Matrix Snm;

    /**
     * Lee el archivo GGM03S.txt y recoge cada fila y lo asigna a Cnm y Snm
     */
void GGM03S();

#endif


