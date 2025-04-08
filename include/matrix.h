#ifndef _MATRIX_
#define _MATRIX_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

    /**
     * @file matrix.h
     * @brief El archivo contiene las funciones
     *  para operaciones con matrices y vectores junto con su definicion
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
class Matrix {
public:
    int n_row, n_column;
	double **data;

    // Parameterized constructor
    /**
     * Crea una matriz [0][n] que simula un vector
     * @param n número de columas del vector
     */
    Matrix(const int n);
    /**
     * Crea una matriz [n_row][n_column]
     * @param n_row número de filas de la matriz
     * @param n_column número de columas de la matriz
     */
    Matrix(const int n_row, const int n_column);
	
	// Member operators
    /**
     * Obtiene el elemento [(n-1)/n_column][(n-1)%n_column]
     * @param n elemnto
     */
	double& operator () (const int n);
    /**
     * Obtiene el elemento [(row-1)][(column-1)]
     * @param n_row fila de la matriz
     * @param n_column columa de la matriz
     */
	double& operator () (const int row, const int column);
    /**
     * Suma 2(this+m) Matrix y devuelve el valor
     * @param m Matrix
     * @return devuelve un Matrix suma de this+ m, sin modificarlos
     */
	Matrix& operator + (Matrix &m);
    /**
     * Resta 2(this-m) Matrix y devuelve el valor
     * @param m Matrix
     * @return devuelve un Matrix resta de this-m, sin modificarlos
     */
	Matrix& operator - (Matrix &m);
    /**
     * Multiplación sobre matrices para  2(this*m) Matrix y devuelve el valor
     * @param m Matrix
     * @return devuelve un Matrix = this*m, sin modificarlos
     */
	Matrix& operator * (Matrix &m);
    /**
     * this*(m^(-1)) y devuelve el valor
     * @param m Matrix
     * @return devuelve un Matrix =this-*m^(-1), sin modificarlos
     */
	Matrix& operator / (Matrix &m);
    /**
     * Se crea una Matrix equivalente a m y se le asigna a this 
     * @param m Matrix
     * @return una Matrix=m, sin modificar las matrices
     */
	Matrix& operator = (Matrix &m);
    /**
     * Suma todas las componentes de this con d y devuelve una nueva Matrix
     * @param d valor a ser operado por cada componente de la matriz
     * @return una Matrix donde todas sus componentes se le suma d, sin modificar las matrices
     */
	Matrix& operator + (double d);
    /**
     * Resta todas las componentes de this con d y devuelve una nueva Matrix
     * @param d valor a ser operado por cada componente de la matriz
     * @return una Matrix donde todas sus componentes se le resta d, sin modificar las matrices
     */
	Matrix& operator - (double d);
    /**
     * Multiplica todas las componentes de this con d y devuelve una nueva Matrix
     * @param d valor a ser operado por cada componente de la matriz
     * @return una Matrix donde todas sus componentes se le multiplica d, sin modificar las matrices
     */
	Matrix& operator * (double d);
    /**
     * Divide todas las componentes de this con d y devuelve una nueva Matrix
     * @param d valor a ser operado por cada componente de la matriz
     * @return una Matrix donde todas sus componentes se le divide d, sin modificar las matrices
     */
	Matrix& operator / (double d);
	
	// Non-member operators
	friend ostream& operator << (ostream &o, Matrix &m);
};

// Operator overloading
ostream& operator << (ostream &o, Matrix &m);


// Methods
    /**
     * Crea una Matrix con todas sus componentes a 0
     * @param n_row numero de filas que tiene la matriz
     * @param n_column numero de columnas que tiene la matriz
     * @return una Matrix tamaño n_row x n_column
     */
	Matrix& zeros(const int n_row, const int n_column);
    /**
     * Crea una Matrix identidad con tamaño size x size
     * @param size dimension de la matriz
     * @return una Matrix tamaño size x size
     */
	Matrix& eye(const int size);
    /**
     * Crea una Matrix traspuesta de m,sin modificar m
     * @param m Matrix
     * @return una Matrix traspuesta de m
     */
	Matrix& transpose(Matrix &m);
    /**
     * Crea una Matrix inversa de m, sin modificar m
     * @param m Matrix que tiene que ser cuadrada, es decir, con el mismo numero de columnas que filas
     * @return una Matrix inversa de m
     */
	Matrix& inv(Matrix &m) ;
    /**
     * Crea una Matrix con todas sus componentes a 0
     * @param n numero de columnas que tiene la matriz
     * @return una Matrix tamaño 1 x n
     */
	Matrix& zeros(const int n);
    /**
     * Devulve la norma 2 de una Matrix que simula un vector
     * @param m Matrix 1 x n_column
     * @return la norma 2 de m
     */
	double norm(Matrix &m) ;
    /**
     * Devulve el producto escalar
     * @param v Matrix con tamaño 1 x 3
     * @param w Matrix con tamaño 1 x 3
     * @return producto escalar de v·w
     */
	double dot(Matrix &v,Matrix &w);
    /**
     * Devulve el producto vectorial
     * @param v Matrix con tamaño 1 x 3
     * @param w Matrix con tamaño 1 x 3
     * @return producto escalar de v x w
     */
	Matrix& cross(Matrix &v,Matrix &w);

    /********************************************************************************************
     * Extrae 
     * @param v Matrix con tamaño 1 x 3
     * @param w Matrix con tamaño 1 x 3
     * @return producto escalar de v x w
     */
    Matrix& extract_vector(Matrix &v,Matrix &w);

    /**
     * Devuelve la matriz resultado de realizar la union entre v y w
     * @param v Matrix con tamaño 1 x n ,n∈N
     * @param w Matrix con tamaño 1 x n ,n∈N
     * @return producto escalar de v x w
     */
    Matrix& union_vector(Matrix &v,Matrix &w);


    /**
     * Extrae la fila i-1 de v y lo devuelve
     * @param v Matrix
     * @param i es la fila tiene que ser >=1 && <=v.n_row
     * @return Matrix con la fila i-1 de v
     */
    Matrix& extract_row(Matrix &v,int i);


    /**
     * Extrae la columna i-1 de v y lo devuelve
     * @param v Matrix
     * @param i es la columna tiene que ser >=1 && <=v.n_column
     * @return Matrix con la columna i-1 de v
     */
    Matrix& extract_column(Matrix &v,int i);


    /**
     * Asigna la fila i-1 con w y lo devuelve
     * @param v Matrix con tamaño x x y ,x∈N y∈N
     * @param w Matrix con tamaño 1 x y ,y∈N
     * @param i es la fila tiene que ser >=1 && <=v.n_row
     * @return Matrix v con la fila i-1 cambiada por los elementos de w
     */
    Matrix& assign_row(Matrix &v,Matrix &w,int i);


    /**
     * Asigna la columna i-1 con w y lo devuelve
     * @param v Matrix con tamaño x x y ,x∈N y∈N
     * @param w Matrix con tamaño 1 x y ,y∈N
     * @param i es la columna tiene que ser >=1 && <=v.n_column
     * @return Matrix v con la columna i-1 cambiada por los elementos de w
     */
    Matrix& assign_column(Matrix &v,Matrix &w,int i);


#endif
