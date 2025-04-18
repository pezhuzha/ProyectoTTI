#include "../include/matrix.h"
#include "../include/R_x.h"
#include "../include/R_y.h"
#include "../include/R_z.h"
#include "../include/AccelPointMass.h"
#include "../include/Cheb3D.h"
#include "../include/EccAnom.h"
#include "../include/Frac.h"
#include "../include/MeanObliquity.h"
#include "../include/Mjday.h"
#include "../include/Mjday_TDB.h"
#include "../include/Position.h"
#include "../include/sign_.h"
#include "../include/timediff.h"
#include "../include/AzElPa.h"
#include "../include/IERS.h"
#include "../include/Legendre.h"
#include "../include/NutAngles.h"
#include "../include/TimeUpdate.h"
#include <cstdio>
#include <cmath>
#include <tuple>

using namespace std;
int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

int m_equals(Matrix A, Matrix B, double p) {
	if (A.n_row != B.n_row || A.n_column != B.n_column)
		return 0;
	else
		for(int i = 1; i <= A.n_row; i++)
			for(int j = 1; j <= A.n_column; j++)
				if(fabs(A(i,j)-B(i,j)) > p) {
					printf("%2.20lf %2.20lf\n",A(i,j),B(i,j));
					return 0;
				}
	
	return 1;
}

int m_sum_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) =  0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = 2; C(1,2) =  2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = 8; C(2,2) = -3; C(2,3) = 1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = -2; C(3,3) = 0; C(3,4) = 7;
	
	Matrix R = A + B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_sub_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) = 0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = -2; C(1,2) = 2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = -6; C(2,2) = 1; C(2,3) = -1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = 4; C(3,3) = 0; C(3,4) = 3;
	
	Matrix R = A - B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_mul_01() {
    int f = 4;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 5; A(1,2) = 8; A(1,3) = 2; A(1,4) = 4;
	A(2,1) = 2; A(2,2) = 2; A(2,3) = 2; A(2,4) = 2;
	A(3,1) = 2; A(3,2) = 2; A(3,3) = 1; A(3,4) = 3;
	A(4,1) = 2; A(4,2) = 2; A(4,3) = 2; A(4,4) = 1;
	
	Matrix B(f, c);

	B(1,1) = 4; B(1,2) = 2; B(1,3) = 9; B(1,4) = 1;
	B(2,1) = 2; B(2,2) = 9; B(2,3) = 2; B(2,4) = 7;
	B(3,1) = 2; B(3,2) = 2; B(3,3) = 4; B(3,4) = 9;
	B(4,1) = 3; B(4,2) = 4; B(4,3) = 2; B(4,4) = 5;
    
	Matrix C(f, c);

	C(1,1) = 52	; C(1,2) = 102	; C(1,3) = 77	; C(1,4) = 99	;
	C(2,1) = 22	; C(2,2) = 34	; C(2,3) = 34	; C(2,4) = 44	;
	C(3,1) = 23	; C(3,2) = 36	; C(3,3) = 32	; C(3,4) = 40	;
	C(4,1) = 19	; C(4,2) = 30	; C(4,3) = 32	; C(4,4) = 39	;

	Matrix R = A * B;
    _assert(m_equals(R, C, 1e-10));
    
    return 0;
}
int m_div_01() {
    int f = 4;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 5; A(1,4) = 5;
	A(2,1) = 2; A(2,2) = 1; A(2,3) = 3; A(2,4) = 6;
	A(3,1) = 5; A(3,2) = 3; A(3,3) = 2; A(3,4) = 3;
	A(4,1) = 1; A(4,2) = 2; A(4,3) = 4; A(4,4) = 1;
	
	Matrix B(f, c);

	B(1,1) = 4; B(1,2) = 2; B(1,3) = 9; B(1,4) = 1;
	B(2,1) = 2; B(2,2) = 9; B(2,3) = 2; B(2,4) = 7;
	B(3,1) = 2; B(3,2) = 2; B(3,3) = 4; B(3,4) = 9;
	B(4,1) = 3; B(4,2) = 4; B(4,3) = 2; B(4,4) = 5;
    
	Matrix C(f, c);

	C(1,1) = 337./963 	; C(1,2) = 350./963 	; C(1,3) = 677./963	; C(1,4) = -271./321;
	C(2,1) = 23./963  	; C(2,2) = -179./963  	; C(2,3) = 592./963	; C(2,4) = 112./321	;
	C(3,1) = 58./963		; C(3,2) = -577./963; C(3,3) = -475./963	; C(3,4) = 743./321	;
	C(4,1) = 430./963 	; C(4,2) = 338./963		; C(4,3) = 98./963	; C(4,4) = -181./321	;

	Matrix R = A / B;
    _assert(m_equals(R, C, 1e-10));
    
    return 0;
}

int m_sum_d_01() {
    int f = 4;
    int c = 4;
	double num=2;
	
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 5; A(1,4) = 5;
	A(2,1) = 2; A(2,2) = 1; A(2,3) = 3; A(2,4) = 6;
	A(3,1) = 5; A(3,2) = 3; A(3,3) = 2; A(3,4) = 3;
	A(4,1) = 1; A(4,2) = 2; A(4,3) = 4; A(4,4) = 1;
    

	Matrix B(f, c);
	B(1,1) = 1+num; B(1,2) = 2+num; B(1,3) = 5+num; B(1,4) = 5+num;
	B(2,1) = 2+num; B(2,2) = 1+num; B(2,3) = 3+num; B(2,4) = 6+num;
	B(3,1) = 5+num; B(3,2) = 3+num; B(3,3) = 2+num; B(3,4) = 3+num;
	B(4,1) = 1+num; B(4,2) = 2+num; B(4,3) = 4+num; B(4,4) = 1+num;


	Matrix R=A+num;

    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int m_sub_d_01() {
    int f = 4;
    int c = 4;
	double num=2;
	
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 5; A(1,4) = 5;
	A(2,1) = 2; A(2,2) = 1; A(2,3) = 3; A(2,4) = 6;
	A(3,1) = 5; A(3,2) = 3; A(3,3) = 2; A(3,4) = 3;
	A(4,1) = 1; A(4,2) = 2; A(4,3) = 4; A(4,4) = 1;


	Matrix B(f, c);
	B(1,1) = 1-num; B(1,2) = 2-num; B(1,3) = 5-num; B(1,4) = 5-num;
	B(2,1) = 2-num; B(2,2) = 1-num; B(2,3) = 3-num; B(2,4) = 6-num;
	B(3,1) = 5-num; B(3,2) = 3-num; B(3,3) = 2-num; B(3,4) = 3-num;
	B(4,1) = 1-num; B(4,2) = 2-num; B(4,3) = 4-num; B(4,4) = 1-num;


	Matrix R=A-num;

    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int m_mul_d_01() {
    int f = 4;
    int c = 4;
	double num=2;

	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 5; A(1,4) = 5;
	A(2,1) = 2; A(2,2) = 1; A(2,3) = 3; A(2,4) = 6;
	A(3,1) = 5; A(3,2) = 3; A(3,3) = 2; A(3,4) = 3;
	A(4,1) = 1; A(4,2) = 2; A(4,3) = 4; A(4,4) = 1;


	Matrix B(f, c);
	B(1,1) = 1.*num; B(1,2) = 2.*num; B(1,3) = 5.*num; B(1,4) = 5.*num;
	B(2,1) = 2.*num; B(2,2) = 1.*num; B(2,3) = 3.*num; B(2,4) = 6.*num;
	B(3,1) = 5.*num; B(3,2) = 3.*num; B(3,3) = 2.*num; B(3,4) = 3.*num;
	B(4,1) = 1.*num; B(4,2) = 2.*num; B(4,3) = 4.*num; B(4,4) = 1.*num;

	Matrix R=A*num;

    _assert(m_equals(R, B, 1e-10));
    
    return 0;
}
int m_div_d_01() {
    int f = 4;
    int c = 4;
	double num=2;

	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 5; A(1,4) = 5;
	A(2,1) = 2; A(2,2) = 1; A(2,3) = 3; A(2,4) = 6;
	A(3,1) = 5; A(3,2) = 3; A(3,3) = 2; A(3,4) = 3;
	A(4,1) = 1; A(4,2) = 2; A(4,3) = 4; A(4,4) = 1;


	Matrix B(f, c);
	B(1,1) = 1./num; B(1,2) = 2./num; B(1,3) = 5./num; B(1,4) = 5./num;
	B(2,1) = 2./num; B(2,2) = 1./num; B(2,3) = 3./num; B(2,4) = 6./num;
	B(3,1) = 5./num; B(3,2) = 3./num; B(3,3) = 2./num; B(3,4) = 3./num;
	B(4,1) = 1./num; B(4,2) = 2./num; B(4,3) = 4./num; B(4,4) = 1./num;

	Matrix R=A/num;

    _assert(m_equals(R, B, 1e-10));
    
    return 0;
}
int m_asig_01() {
    int f = 4;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 5; A(1,4) = 5;
	A(2,1) = 2; A(2,2) = 1; A(2,3) = 3; A(2,4) = 6;
	A(3,1) = 5; A(3,2) = 3; A(3,3) = 2; A(3,4) = 3;
	A(4,1) = 1; A(4,2) = 2; A(4,3) = 4; A(4,4) = 1;
    
	Matrix B=A;

    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}
int m_zeros_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 0; A(1,3) = 0; A(1,4) = 0;
	A(2,1) = 0; A(2,2) = 0; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 0; A(3,4) = 0;
	
	Matrix B = zeros(3, 4);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_eye_01() {
    int f = 3;
	
	
	Matrix A(f, f);
	A(1,1) = 1; A(1,2) = 0; A(1,3) = 0;
	A(2,1) = 0; A(2,2) = 1; A(2,3) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 1;

	
	Matrix B = eye(f);

    _assert(m_equals(A, B, 1e-10));
    
    
    return 0;
}

int m_transpose_01() {
    int f = 3;
    int c = 3;
	
	
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 4; A(1,3) = 9;
	A(2,1) = 2; A(2,2) = 3; A(2,3) = 8;
	A(3,1) = 5; A(3,2) = 6; A(3,3) = 7;

	
	Matrix B(f,c);

	B(1,1) = 1; B(1,2) = 2; B(1,3) = 5;
	B(2,1) = 4; B(2,2) = 3; B(2,3) = 6;
	B(3,1) = 9; B(3,2) = 8; B(3,3) = 7;

   	Matrix R=transpose(A);

    _assert(m_equals(R, B, 1e-10));
    
    return 0;
}
int m_inv_01() {
    int f = 2;
	
	
	Matrix A(f, f);
	A(1,1) = 5; A(1,2) = 2;
	A(2,1) =-7; A(2,2) = -3;

	
	Matrix B(f,f);

	B(1,1) = 3; B(1,2) = 2;
	B(2,1) = -7; B(2,2) = -5;
   	Matrix R=inv(A);
    _assert(m_equals(R, B, 1e-10));
    return 0;
}
int m_inv_02() {
    int f = 4;
	
	
	Matrix A(f, f);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 5; A(1,4) = 5;
	A(2,1) = 2; A(2,2) = 1; A(2,3) = 3; A(2,4) = 6;
	A(3,1) = 5; A(3,2) = 3; A(3,3) = 2; A(3,4) = 3;
	A(4,1) = 1; A(4,2) = 2; A(4,3) = 4; A(4,4) = 1;
	
	Matrix B(f, f);

	B(1,1)=-47./42; B(1,2) = 5./6; B(1,3) = -1./14 ; B(1,4) =17./21;
	B(2,1)=41./21; B(2,2) = -5./3; B(2,3) = 4./7 ; B(2,4) =-31./21;
	B(3,1) = -17./21 ; B(3,2) = 2./3 ; B(3,3) = -2./7; B(3,4) =19./21;
	B(4,1) = 19./42 ; B(4,2) = -1./6; B(4,3) = 1./14 ; B(4,4) = -10./21;

   	Matrix R=inv(A);
    _assert(m_equals(R, B, 1e-10));
    return 0;
}
int m_norm_01() {
    int f = 4;
	
	
	Matrix A(f);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 5; A(1,4) = 5;
	
	double B=7.4161984870956629487113974408007;

   	double R=norm(A);
    _assert(fabs(R-B)<1e-10);
    return 0;
}

int m_dot_01() {
	int f = 3;
	Matrix A(f);
	A(1,1) = 2; A(1,2) = 1; A(1,3) = 0;
	
	Matrix B(f);
	B(1,1)= 3; B(1,2) = 5; B(1,3) = 6 ;
	
	
	double R=11;
	
	double C=dot(A,B);
	
    _assert(fabs(R-C)<1e-10);
    
    return 0;
}
int m_cross_01() {
	int f = 3;
	Matrix A(f);
	Matrix B(f);
	A(1,1) = 2; A(1,2) = 1; A(1,3) = 0;
	B(1,1)= 3; B(1,2) = 5; B(1,3) = 6 ;
    
	Matrix R(f);
	R(1,1)= 6; R(1,2) = -12; R(1,3) = 7 ;
	
	Matrix C=cross(A,B);
	
    _assert(m_equals(R, C, 1e-10));
    
    return 0;
}
int m_extract_vector_01() {
	int f = 5;
	Matrix A(f);
	
	A(1,1) = 2; A(1,2) = 1; A(1,3) = 0; A(1,4) = 5; A(1,5) = 1;

	Matrix B=extract_vector(A,1,3);
    
	Matrix R(3);
	R(1,1)= 2; R(1,2) = 1; R(1,3) = 0 ;

	
    _assert(m_equals(R, B, 1e-10));
    
    return 0;
}

int m_union_vector_01() {
	int f = 3;
	Matrix A(f);
	Matrix B(f);
	A(1,1) = 2; A(1,2) = 1; A(1,3) = 0;
	B(1,1)= 3; B(1,2) = 1; B(1,3) = 6 ;
    
	Matrix R(5);
	R(1,1) = 2; R(1,2) = 1; R(1,3) = 0;
	R(1,4)= 3; R(1,5) = 6 ;
	
	Matrix C=union_vector(A,B);
    _assert(m_equals(R, C, 1e-10));
    
    return 0;
}

int m_extract_row_01() {
	int f = 3;
	Matrix A(f,f);
	
	A(1,1) = 2; A(1,2) = 1; A(1,3) = 0;
	A(2,1) = 2; A(2,2) = 1; A(2,3) = 3;
	A(3,1) = 5; A(3,2) = 3; A(3,3) = 2;

	Matrix B=extract_row(A,f);
    
	Matrix R(f);
	R(1,1)= 5; R(1,2) = 3; R(1,3) = 2 ;

    _assert(m_equals(R, B, 1e-10));
    
    return 0;

}

int m_extract_column_01() {
	int f = 3;
	Matrix A(f,f);
	
	A(1,1) = 2; A(1,2) = 1; A(1,3) = 0;
	A(2,1) = 2; A(2,2) = 1; A(2,3) = 3;
	A(3,1) = 5; A(3,2) = 3; A(3,3) = 2;

	Matrix B=extract_column(A,f);
    
	Matrix R(f);
	R(1,1)= 0; R(1,2) = 3; R(1,3) = 2 ;

	
    _assert(m_equals(R, B, 1e-10));
    
    return 0;
}

int m_assign_row_01() {
	int f = 3;
	Matrix A(f,f);
	
	A(1,1) = 2; A(1,2) = 1; A(1,3) = 0;
	A(2,1) = 2; A(2,2) = 1; A(2,3) = 3;
	A(3,1) = 5; A(3,2) = 3; A(3,3) = 2;

	Matrix B(f);

	B(1,1)= 3; B(1,2) = 5; B(1,3) = 6 ;


	Matrix C=assign_row(A,B,f);
    
	Matrix R(f,f);
	R(1,1) = 2; R(1,2) = 1; R(1,3) = 0;
	R(2,1) = 2; R(2,2) = 1; R(2,3) = 3;
	R(3,1) = 3; R(3,2) = 5; R(3,3) = 6;

    _assert(m_equals(R, C, 1e-10));
    
    return 0;
}

int m_assign_column_01() {
	int f = 3;
	Matrix A(f,f);
	
	A(1,1) = 2; A(1,2) = 1; A(1,3) = 0;
	A(2,1) = 2; A(2,2) = 1; A(2,3) = 3;
	A(3,1) = 5; A(3,2) = 3; A(3,3) = 2;

	Matrix B(f);

	B(1,1)= 3; B(1,2) = 5; B(1,3) = 6 ;

	Matrix C=assign_column(A,B,f);

	Matrix R(f,f);
	R(1,1) = 2; R(1,2) = 1; R(1,3) = 3;
	R(2,1) = 2; R(2,2) = 1; R(2,3) = 5;
	R(3,1) = 5; R(3,2) = 3; R(3,3) = 6;

	
    _assert(m_equals(R, C, 1e-10));
    
    return 0;
}

int m_R_x_01() {	
	Matrix A=R_x(10);


	Matrix R(3,3);                 
	R(1,1) = 1; R(1,2) = 0; R(1,3) = 0;
	R(2,1) = 0; R(2,2) = -0.839071529076452; R(2,3) = -0.54402111088937;
	R(3,1) = 0; R(3,2) = 0.54402111088937; R(3,3) = -0.839071529076452;

    _assert(m_equals(R, A, 1e-10));
    
    return 0;
}

int m_R_y_01() {
	
	Matrix A=R_y(10);

	Matrix R(3,3);                
	R(1,1) = -0.839071529076452 ; R(1,2) = 0; R(1,3) = 0.54402111088937;
	R(2,1) = 0					; R(2,2) = 1; R(2,3) = 0;
	R(3,1) = -0.54402111088937  ; R(3,2) = 0; R(3,3) = -0.839071529076452;

	
    _assert(m_equals(R, A, 1e-10));
    
    return 0;
}

int m_R_z_01() {
	
	Matrix A=R_z(10);


	Matrix R(3,3);                   
	R(1,1) = -0.839071529076452; R(1,2) = -0.54402111088937; R(1,3) = 0;
	R(2,1) = 0.54402111088937; R(2,2) = -0.839071529076452; R(2,3) = 0;
	R(3,1) = 0; R(3,2) = 0; R(3,3) = 1;


    _assert(m_equals(R, A, 1e-10));
    
    return 0;
}

int m_AccelPointMass_01() {
	
	Matrix A(3);    
	
	A(1,1) = 1; A(1,2) = 1; A(1,3) = 1;

	Matrix B(3);

	B(1,1)= 2; B(1,2) = 3; B(1,3) = 4 ;

	Matrix C=AccelPointMass(A,B,10);
	
	Matrix R(3);  
	R(1,1) = 0.0628351366133708; R(1,2) = 0.189703148460208; R(1,3) = 0.316571160307045;
	
    _assert(m_equals(R, A, 1e-10));
    
    return 0;
}
int m_Cheb3D_01() {
	double f = 3;
	
	Matrix A(f);
	
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;

	Matrix B(f);

	B(1,1)= 1; B(1,2) = 2; B(1,3) = 3;
	
	Matrix C(f);

	C(1,1)= 5; C(1,2) = 2; C(1,3) = 3;

	Matrix D=Cheb3D(1,3,0.5,1,A,B,C);
	
	Matrix R(f);
	R(1,1) = 6; R(1,2) = 6; R(1,3) = 10;
	
    _assert(m_equals(R, D, 1e-10));
    
    return 0;
}
int m_EccAnom_01() {
	
	double R = 2.38006127313934;
	double D=EccAnom(1,2);
	
    _assert(fabs(R-D)< 1e-10);
    
    return 0;
}
int m_Frac_01() {
	
	double R = 0.3801;
	double D=Frac(2.3801);
	
    _assert(fabs(R-D)< 1e-10);
    
    return 0;
}

int m_MeanObliquity_01() {
	
	double R = 0.409412815476201;
	double D=MeanObliquity(41);
	
    _assert(fabs(R-D)< 1e-10);
    
    return 0;
}
//********************************************************************************************************************************************************************************************************

int m_Mjday_01() {
	
	double R = 0.409412815476201;
	double D= Mjday(2025,5,5);
	
    _assert(fabs(R-D)< 1e-10);
    
    return 0;
}
int m_Mjday_TDB_01() {
	
	double R = 0.409412815476201;
	double D= Mjday_TDB(2025);
	
    _assert(fabs(R-D)< 1e-10);
    
    return 0;
}
int m_Position_01() {
	
	Matrix R(3);
	R(1) = 6; R(2) = 6; R(3) = 10;

	Matrix D= Position(2,3,4);
	
    _assert(m_equals(R, D, 1e-10));
    
    return 0;
}
int m_sign__01() {
	
	double R = -4;
	double D= sign_(4,-3);
	
    _assert(fabs(R-D)< 1e-10);
    
    return 0;
}

int m_timediff_01() {
	
	double R0 = -4;
	double R1 = -4;
	double R2 = -4;
	double R3 = -4;
	double R4 = -4;
	auto D= timediff(4,10);
	_assert(fabs(get<0>(D)-R0)< 1e-10);
	_assert(fabs(get<1>(D)-R1)< 1e-10);
	_assert(fabs(get<2>(D)-R2)< 1e-10);
	_assert(fabs(get<3>(D)-R3)< 1e-10);
	_assert(fabs(get<4>(D)-R4)< 1e-10);
    
    return 0;
}
int m_AzElPa_01() {
	Matrix A(3);
	
	A(1) = 1; A(2) = 2; A(3) = 3;
	
	double R0=0;

	double R1=1;

	Matrix R2(3);
	
	R2(1) = 1; R2(2) = 2; R2(3) = 3;
	Matrix R3(3);
	
	R3(1) = 1; R3(2) = 2; R3(3) = 3;

	auto D= AzElPa(A);
	_assert(fabs(get<0>(D)-R0)< 1e-10);
	_assert(fabs(get<1>(D)-R1)< 1e-10);
	_assert(m_equals(get<2>(D),R2,1e-10));
	_assert(m_equals(get<3>(D),R3,1e-10));
	
    
    return 0;
}
int m_IERS_01() {

	
	Matrix A(4, 13);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 5; A(1,4) = 5;
	A(2,1) = 2; A(2,2) = 1; A(2,3) = 3; A(2,4) = 6;
	A(3,1) = 5; A(3,2) = 3; A(3,3) = 2; A(3,4) = 3;
	A(4,1) = 1; A(4,2) = 2; A(4,3) = 4; A(4,4) = 1;
	A(1,5) = 1; A(1,7) = 2; 
	A(2,5) = 2; A(2,7) = 1; 
	A(3,5) = 5; A(3,7) = 3; 
	A(4,5) = 1; A(4,7) = 2; 
	A(1,6) = 1; A(1,8) = 2; 
	A(2,6) = 2; A(2,8) = 1; 
	A(3,6) = 5; A(3,8) = 3;
	A(4,6) = 1; A(4,8) = 2; 
	A(1,9) = 5; 
	A(2,9) = 3; 
	A(3,9) = 2;
	A(4,9) = 4; 
	A(1,10)= 5; 
	A(2,10)= 3; 
 	A(3,10)= 2; 
	A(4,10)= 4;
	A(1,11) = 5;
	A(2,11) = 6; 
	A(3,11) = 3;
	A(4,11) = 1;
	A(1,12) = 5;
	A(2,12) = 6;
 	A(3,12) = 3;
	A(4,12) = 1;
	A(1,13)= 5; 
	A(2,13)= 3; 
 	A(3,13)= 2;
	A(4,13)= 4;

	double R0 = -4;
	double R1 = -4;
	double R2 = -4;
	double R3 = -4;
	double R4 = -4;
	double R5 = -4;
	double R6 = -4;
	double R7 = -4;
	double R8 = -4;


	auto D= IERS(A,1,'l');
	_assert(fabs(get<0>(D)-R0)< 1e-10);
	_assert(fabs(get<1>(D)-R1)< 1e-10);
	_assert(fabs(get<2>(D)-R2)< 1e-10);
	_assert(fabs(get<3>(D)-R3)< 1e-10);
	_assert(fabs(get<4>(D)-R4)< 1e-10);
	_assert(fabs(get<5>(D)-R5)< 1e-10);
	_assert(fabs(get<6>(D)-R6)< 1e-10);
	_assert(fabs(get<7>(D)-R7)< 1e-10);
	_assert(fabs(get<8>(D)-R8)< 1e-10);
	
    
    return 0;
}

int m_Legendre_01() {


	Matrix R0 (3,3);                
	R0(1,1) = -0.839071529076452 ; R0(1,2) = 0; R0(1,3) = 0.54402111088937;
	R0(2,1) = 0					; R0(2,2) = 1; R0(2,3) = 0;
	R0(3,1) = -0.54402111088937  ; R0(3,2) = 0; R0(3,3) = -0.839071529076452;
	Matrix R1 (3,3);                
	R1(1,1) = -0.839071529076452 ; R1(1,2) = 0; R1(1,3) = 0.54402111088937;
	R1(2,1) = 0					; R1(2,2) = 1; R1(2,3) = 0;
	R1(3,1) = -0.54402111088937  ; R1(3,2) = 0; R1(3,3) = -0.839071529076452;

	auto D= Legendre(3,3,1);
	_assert(m_equals(get<0>(D),R0, 1e-10));
	_assert(m_equals(get<1>(D),R1,1e-10));
	
    
    return 0;
}
int m_NutAngles_01() {


	double R0 = -4;
	double R1 = -4;

	auto D= NutAngles(3);

	_assert(fabs(get<0>(D)-R0)< 1e-10);
	_assert(fabs(get<1>(D)-R1)< 1e-10);
	
    
    return 0;
}
int m_TimeUpdate_01() {

	Matrix A(3, 3);
	A(1,1) = 1; A(1,2) = 4; A(1,3) = 9;
	A(2,1) = 2; A(2,2) = 3; A(2,3) = 8;
	A(3,1) = 5; A(3,2) = 6; A(3,3) = 7;

	
	Matrix B(3,3);

	B(1,1) = 1; B(1,2) = 2; B(1,3) = 5;
	B(2,1) = 4; B(2,2) = 3; B(2,3) = 6;
	B(3,1) = 9; B(3,2) = 8; B(3,3) = 7;

	Matrix C(3, 3);

	C(1,1) = 52	; C(1,2) = 102	; C(1,3) = 77	;
	C(2,1) = 22	; C(2,2) = 34	; C(2,3) = 34	;
	C(3,1) = 23	; C(3,2) = 36	; C(3,3) = 32	;


	Matrix D(3, 3);

	D(1,1) = 52	; D(1,2) = 102	; D(1,3) = 77	;
	D(2,1) = 22	; D(2,2) = 34	; D(2,3) = 34	;
	D(3,1) = 23	; D(3,2) = 36	; D(3,3) = 32	;

	Matrix R = TimeUpdate(A,B,C);
	
	_assert(m_equals(R,D,1e-10));
    
    return 0;
}
//********************************************************************************************************************************************************************************************************
int all_tests()
{
    _verify(m_sum_01);
    _verify(m_sub_01);
    _verify(m_mul_01);
    _verify(m_div_01);
    _verify(m_sum_d_01);
    _verify(m_sub_d_01);
    _verify(m_mul_d_01);
    _verify(m_div_d_01);
    _verify(m_asig_01);
    _verify(m_zeros_01);
    _verify(m_eye_01);
    _verify(m_transpose_01);
    _verify(m_inv_01);
    _verify(m_inv_02);
    _verify(m_norm_01);
    _verify(m_dot_01);
    _verify(m_cross_01);
    _verify(m_extract_vector_01);
    _verify(m_union_vector_01);
    _verify(m_extract_row_01);
    _verify(m_extract_column_01);
    _verify(m_assign_row_01);
    _verify(m_assign_column_01);
    _verify(m_R_x_01);
    _verify(m_R_y_01);
    _verify(m_R_z_01);
    _verify(m_Cheb3D_01);
    _verify(m_EccAnom_01);
    _verify(m_Frac_01);
    _verify(m_MeanObliquity_01);
    _verify(m_Mjday_01);
    _verify(m_Mjday_TDB_01);
    _verify(m_Position_01);
    _verify(m_sign__01);
    _verify(m_timediff_01);
    _verify(m_AzElPa_01);
    _verify(m_IERS_01);
    _verify(m_Legendre_01);
    _verify(m_NutAngles_01);


    return 0;
}


int main()
{
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return (result != 0);
}
