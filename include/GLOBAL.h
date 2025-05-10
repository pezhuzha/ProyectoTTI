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

typedef struct{
	double Mjd_UTC,Mjd_TT;
	int n,m,sun,moon,planets;
} Param;

extern Param AuxParam;
extern Matrix eopdata;
extern Matrix Cnm;
extern Matrix Snm;
extern Matrix PC;

    /**
     * Lee el archivo eop19620101.txt y recoge cada fila y lo asigna a eopdata
     * @param c    número de filas a recoger
     */
void eop19620101(int c=21413);

    /**
     * Lee el archivo GGM03S.txt y recoge cada fila y lo asigna a Cnm y Snm
     * @param c    dimension de la matriz
     */
void GGM03S(int n=181);

    /**
     * Lee el archivo DE430Coeff.txt y recoge cada fila y lo asigna a PC
     * @param row    número de filas a recoger
     * @param column    número de columnas a recoger
     */
void DE430Coeff(int row=2285,int column=1020);

    /**
     * Carga AuxParam
     */
void AuxParamLoad();

#endif


