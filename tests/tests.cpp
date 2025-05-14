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
#include "../include/GLOBAL.h"
#include "../include/AccelHarmonic.h"
#include "../include/EqnEquinox.h"
#include "../include/JPL_Eph_DE430.h"
#include "../include/LTC.h"
#include "../include/NutMatrix.h"
#include "../include/PoleMatrix.h"
#include "../include/PrecMatrix.h"
#include "../include/gmst.h"
#include "../include/gast.h"
#include "../include/MeasUpdate.h"
#include "../include/G_AccelHarmonic.h"
#include "../include/GHAMatrix.h"
#include "../include/Accel.h"
#include "../include/VarEqn.h"
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
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
	
	Matrix B(f);
	B(1,1)= 1; B(1,2) = 2; B(1,3) = 3 ;
	
	
	double R=14;
	
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
    
	Matrix R(6);
	R(1,1) = 2; R(1,2) = 1; R(1,3) = 0;
	R(1,4)= 3; R(1,5) = 1;R(1,6) = 6 ;
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


int m_Mjday_01() {
	
	double R = 60800;
	double D= Mjday(2025,5,5);
	
    _assert(fabs(R-D)< 1e-10);
    
    return 0;
}
int m_Mjday_TDB_01() {
	
	double R = 2025.0000000092;
	double D= Mjday_TDB(2025);
	
    _assert(fabs(R-D)< 1e-10);
    
    return 0;
}
int m_Position_01() {
	
	Matrix R(3);
	R(1) = 2.627855739427486e+06; R(2) = -5.741969545549633e+06; R(3) = 8.941173180321892e+05;
	
	Matrix D= Position(2,3,4);
    _assert(m_equals(R, D, 1e-8));
    
    return 0;
}
int m_sign__01() {
	
	double R = -4;
	double D= sign_(4,-3);
	
    _assert(fabs(R-D)< 1e-10);
    
    return 0;
}
int m_timediff_01() {
	
	double R0 = -6;
	double R1 = 9;
	double R2 = 13;
	double R3 = 42.184;
	double R4 = -9;
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
	
	double R0=0.463647609000806;

	double R1=0.930274014115472;

	Matrix R2(3);
	
	R2(1) = 0.4; R2(2) = -0.2; R2(3) = 0;
	Matrix R3(3);
	
	R3(1) = -0.095831484749991; R3(2) = -0.191662969499982; R3(3) = 0.159719141249985;

	auto [Az, El, dAds, dEds]= AzElPa(A);
	_assert(fabs(Az-R0)< 1e-10);
	_assert(fabs(El-R1)< 1e-10);
	_assert(m_equals(dAds,R2,1e-10));
	_assert(m_equals(dEds,R3,1e-10));
	
    
    return 0;
}

int m_IERS_01() {

	
	eop19620101(21413);

	double R0 = -5.59518621231704e-07;
	double R1 = 2.33458634442529e-06;
	double R2 =  0.3260677;
	double R3 = 0.0027213;
	double R4 = -1.16864337831454e-07;
	double R5 = -2.48709418409192e-08;
	double R6 = -8.19335121075116e-10;
	double R7 =  -1.53201123230613e-09;
	double R8 = 29;


	auto [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC]= IERS(eopdata,49746,'l');
	_assert(fabs(x_pole-R0)< 1e-10);
	_assert(fabs(y_pole-R1)< 1e-10);
	_assert(fabs(UT1_UTC-R2)< 1e-10);
	_assert(fabs(LOD-R3)< 1e-10);
	_assert(fabs(dpsi-R4)< 1e-10);
	_assert(fabs(deps-R5)< 1e-10);
	_assert(fabs(dx_pole-R6)< 1e-10);
	_assert(fabs(dy_pole-R7)< 1e-10);
	_assert(fabs(TAI_UTC-R8)< 1e-10);
	
    
    return 0;
}

int m_Legendre_01() {


	Matrix R0 (4,4);                
	R0(1,1) = 1; R0(1,2) = 0; R0(1,3) = 0; R0(1,4) = 0;
	R0(2,1) =  1.4574704987823; R0(2,2) =  0.935831045210238; R0(2,3) = 0; R0(2,4) = 0;
	R0(3,1) =  1.25691645573063 ; R0(3,2) = 1.76084689542256  ; R0(3,3) = 0.565313394670859 ; R0(3,4) = 0;
	R0(4,1) = 0.601515831515714; R0(4,2) = 2.22381140389174 ; R0(4,3) = 1.25857019087392 ; R0(4,4) = 0.329913047636197;
	Matrix R1 (4,4);                
	R1(1,1) = 0; R1(1,2) = 0; R1(1,3) = 0; R1(1,4) = 0;
	R1(2,1) = 0.935831045210238 ; R1(2,2) = -1.4574704987823 ; R1(2,3) = 0; R1(2,4) = 0;
	R1(3,1) =  3.0498762872218; R1(3,2) = -1.61172976752398 ; R1(3,3) = -1.76084689542256  ; R1(3,4) = 0;
	R1(4,1) = 5.44720322371707  ; R1(4,2) =  0.516567339757783 ; R1(4,3) =  -3.11209524837966 ; R1(4,4) = -1.54142738655916;

	auto [pnm, dpnm]= Legendre(3,3,1);
	_assert(m_equals(pnm,R0, 1e-10));
	_assert(m_equals(dpnm,R1,1e-10));
	
    
    return 0;
}
int m_NutAngles_01() {


	double R0 = 2.72256565175042e-05;
	double R1 =  3.87947551912632e-05;

	auto [dpsi, deps]= NutAngles(3);
	_assert(fabs(dpsi-R0)< 1e-10);
	_assert(fabs(deps-R1)< 1e-10);
	
    
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

	D(1,1) = 462	; D(1,2) = 702	; D(1,3) = 1087	;
	D(2,1) = 694	; D(2,2) = 989	; D(2,3) = 1596	;
	D(3,1) = 1257	; D(3,2) = 1746	; D(3,3) = 2746	;

	Matrix R = TimeUpdate(A,B,C);
	
	_assert(m_equals(R,D,1e-10));
    
    return 0;
}
int m_AccelHarmonic_01() {

	Matrix R0(3);
	R0(1)=2.42488766379455e+34;
	R0(2)=2.65762182552943e+34;
	R0(3)=2.65762182552943e+34;
	
	Matrix A(3);
	A(1)=1.0;
	A(2)=2.0;
	A(3)=3.0;
	A=transpose(A);
	Matrix B(3,3);
	B(1,1)= 1.0; B(1,2) = 1.0; B(1,3) = 1.0 ;
	B(2,1)= 1.0; B(2,2) = 2.0; B(2,3) = 2.0 ;
	B(3,1)= 1.0; B(3,2) = 3.0; B(3,3) = 3.0 ;
	
	Matrix R = AccelHarmonic(A,B,5,5);
	_assert(m_equals(R,R0,R0(1)*1e-10));
    
    return 0;
}

int m_EqnEquinox_01() {

	double R0=2.6045897022442e-05;
	double R = EqnEquinox(5);
	
	_assert(fabs(R-R0)< 1e-10);
    
    return 0;
}

int m_JPL_Eph_DE430_01() {

	Matrix R0(3);
	R0(1)=147208460159.245;
	R0(2)=  54592844683.9181;
	R0(3)= 15319523517.8098;
	Matrix R1(3);
	R1(1)= 72752904522.483;
	R1(2)=2340227175.73022;
	R1(3)=1670913926.26657;
	Matrix R2(3);
	R2(1)=-108493583087.765;
	R2(2)=-97599455066.7732;
	R2(3)=-42280555048.3341;
	Matrix R3(3);
	R3(1)=-131548434829.954;
	R3(2)= 156673495960.063;
	R3(3)= 75876841465.0958;
	Matrix R4(3);
	R4(1)=125152801267.633;
	R4(2)= 801496792415.556;
	R4(3)=343590087271.679;
	Matrix R5(3);
	R5(1)= 1532648188929.54;
	R5(2)=-30170201818.4012;
	R5(3)=-71835402852.1091;
	Matrix R6(3);
	R6(1)= 1707620424797.68;
	R6(2)= 2344787198692.63;
	R6(3)=1003869836966.23;
	Matrix R7(3);
	R7(1)= 4578089056071.27;
	R7(2)= 104288160816.022;
	R7(3)= -66258775563.9107;
	Matrix R8(3);
	R8(1)=2885822907234;
	R8(2)= -3876433792795.09;
	R8(3)=-2034695099132.45;
	Matrix R9(3);
	R9(1)=-295931131.772483;
	R9(2)= 224622149.622798;
	R9(3)= 120216992.495936;
	Matrix R10(3);
	R10(1)= 107770931597.498;
	R10(2)=96865992179.6201;
	R10(3)=41989334136.8052;
	auto [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,r_Sun] = JPL_Eph_DE430(60800);
	_assert(m_equals(r_Mercury,R0, abs(R0(1)*1e-10)));
	_assert(m_equals(r_Venus,R1,abs(R1(1)*1e-10)));
	_assert(m_equals(r_Earth,R2, abs(R2(1)*1e-10)));
	_assert(m_equals(r_Mars,R3,abs(R3(1)*1e-10)));
	_assert(m_equals(r_Jupiter,R4, abs(R4(1)*1e-10)));
	_assert(m_equals(r_Saturn,R5,abs(R5(1)*1e-10)));
	_assert(m_equals(r_Uranus,R6, abs(R6(1)*1e-10)));
	_assert(m_equals(r_Neptune,R7,abs(R7(1)*1e-10)));
	_assert(m_equals(r_Pluto,R8, abs(R8(1)*1e-10)));
	_assert(m_equals(r_Moon,R9,abs(R9(1)*1e-10)));
	_assert(m_equals(r_Sun,R10, abs(R10(1)*1e-10)));
    
    return 0;
}


int m_LTC_01() {

	Matrix A=LTC(10,10);


	Matrix R(3,3);                   
	R(1,1) = 0.54402111088937; R(1,2) = -0.839071529076452; R(1,3) = 0;
	R(2,1) = -0.456472625363814; R(2,2) = -0.295958969093304; R(2,3) = -0.839071529076452;
	R(3,1) = 0.704041030906696; R(3,2) = 0.456472625363814; R(3,3) = -0.54402111088937;


    _assert(m_equals(R, A, 1e-10));
    
    return 0;
}


int m_NutMatrix_01() {

	Matrix A=NutMatrix(10);


	Matrix R(3,3);                   
	R(1,1) = 0.999999999492159; R(1,2) = -2.92358797494727e-05 ; R(1,3) = -1.26864277798066e-05;
	R(2,1) = 2.9235393806329e-05 ; R(2,2) = 0.999999998839099; R(2,3) = -3.83026695232602e-05;
	R(3,1) = 1.26875475773192e-05; R(3,2) =  3.83022986110704e-05 ; R(3,3) =  0.99999999918598;


    _assert(m_equals(R, A, 1e-10));
    
    return 0;
}
int m_PoleMatrix_01() {

	Matrix A=PoleMatrix(10,10);


	Matrix R(3,3);                   
	R(1,1) = -0.839071529076452; R(1,2) = 0.295958969093304; R(1,3) = 0.456472625363814;
	R(2,1) = 0; R(2,2) = -0.839071529076452; R(2,3) = 0.54402111088937;
	R(3,1) = 0.54402111088937; R(3,2) = 0.456472625363814; R(3,3) = 0.704041030906696;


    _assert(m_equals(R, A, 1e-10));
    
    return 0;
   }

int m_PrecMatrix_01() {

	Matrix A=PrecMatrix(100,1);


	Matrix R(3,3);                   
	R(1,1) = 0.999999997819034; R(1,2) =  6.05590736738844e-05; R(1,3) = 2.63539319986234e-05;
	R(2,1) =  -6.05590736738844e-05; R(2,2) = 0.999999998166299 ; R(2,3) = -7.97984483806257e-10;
	R(3,1) = -2.63539319986234e-05; R(3,2) =   -7.97985227435292e-10; R(3,3) =  0.999999999652735;


    _assert(m_equals(R, A, 1e-10));
    
    return 0;
   }
   int m_gmst_01() {

	double A=gmst(10);


	double R=1.14523606099042;


    _assert(fabs(R-A)< 1e-10);
    
    return 0;
   }

   int m_gast_01() {

	double A=gast(10);


	double R=1.14526529687017;


    _assert(fabs(R-A)< 1e-10);
    
    return 0;
   }

int m_MeasUpdate_01() {

	Matrix A(3);
	A(1)=1;
	A(2)=2;
	A(3)=3;
	Matrix B=transpose(A);

	Matrix C(3,3);                   
	C(1,1) = 1; C(1,2) = 2; C(1,3) = 3;
	C(2,1) = 6; C(2,2) = 2; C(2,3) = 3;
	C(3,1) = 8; C(3,2) = 2; C(3,3) = 3;

	auto [K, x, P]=MeasUpdate(B,2,3,4,A,C,3);


	Matrix R0(3);
	R0(1)=0.106870229007634;
	R0(2)=  0.145038167938931;
	R0(3)= 0.16030534351145;
	R0=transpose(R0);
	Matrix R1(3);
	R1(1)= 0.893129770992366;
	R1(2)= 1.85496183206107;
	R1(3)=2.83969465648855;
	R1=transpose(R1);

	Matrix R2(3,3);                   
	R2(1,1) = -2.95419847328244  ; R2(1,2) = 0.717557251908397; R2(1,3) = 1.0763358778626;
	R2(2,1) = 0.633587786259541; R2(2,2) = 0.259541984732824; R2(2,3) = 0.389312977099237;
	R2(3,1) = 2.06870229007634; R2(3,2) = 0.0763358778625954 ; R2(3,3) = 0.114503816793893;
    _assert(m_equals(R0, K, 1e-10));
    _assert(m_equals(R1, x, 1e-10));
    _assert(m_equals(R2, P, 1e-10));
    
    return 0;
   }
   int m_G_AccelHarmonic_01() {

	Matrix R2(3,3);                   
	R2(1,1) = -2.0122905124052e+34   ; R2(1,2) =-3.29511632095482e+34 ; R2(1,3) = -3.29511632095482e+34;
	R2(2,1) = -3.11072371483497e+34  ; R2(2,2) = -3.38401367341354e+34 ; R2(2,3) = -3.38401367341354e+34;
	R2(3,1) =  -3.11072371483497e+34 ; R2(3,2) =  -3.38401367341354e+34 ; R2(3,3) = -3.38401367341354e+34;
	
	Matrix A(3);
	A(1)=1.0;
	A(2)=2.0;
	A(3)=3.0;
	A=transpose(A);
	Matrix B(3,3);
	B(1,1)= 1.0; B(1,2) = 1.0; B(1,3) = 1.0 ;
	B(2,1)= 1.0; B(2,2) = 2.0; B(2,3) = 2.0 ;
	B(3,1)= 1.0; B(3,2) = 3.0; B(3,3) = 3.0 ;
	
	Matrix R = G_AccelHarmonic(A,B,5,5);
	_assert(m_equals(R,R2,fabs(R2(1)*1e-10)));
    
    return 0;
}
int m_GHAMatrix_01() {

	Matrix R(3,3);                   
	R(1,1) = 0.412804512414729; R(1,2) = 0.910819649837463 ; R(1,3) = 0;
	R(2,1) = -0.910819649837463 ; R(2,2) = 0.412804512414729; R(2,3) = 0;
	R(3,1) = 0; R(3,2) = 0; R(3,3) = 1;
	
	
	Matrix A = GHAMatrix(10);

	_assert(m_equals(R,A,1e-10));
    
    return 0;
}

int m_Accel_01() {

	Matrix R(6);
	R(1)=1.0;
	R(2)=2.0;
	R(3)=3.0;
	R(4)=-9.52489066332755e+131;
	R(5)=-1.68703107956274e+132;
	R(6)=-4.07471909292663e+132;
	
	Matrix A(6);
	A(1)=1.0;
	A(2)=2.0;
	A(3)=3.0;
	A(4)=1.0;
	A(5)=2.0;
	A(6)=3.0;
	A=transpose(A);
	Matrix B = Accel(10,A);


	_assert(m_equals(R,B,abs(R(6)*1e-10)));
    
    return 0;
}
int m_VarEqn_01() {
Matrix A(42);
	A(1)=7101800.90695315;
	A(2)=1293997.58115302;
	A(3)=10114.014948955;
	A(4)= 573.068082065557;
	A(5)= -3085.15736953138;
	A(6)=      -6736.03068347156;
	A(7)=        1.0000293469741;
	A(8)=    8.22733917593032e-06;
	A(9)=   2.17104932968693e-07;
	A(10)=    1.08925458231315e-05;
	A(11)=   3.04673932160225e-06;
	A(12)=    6.63504292706821e-08;
	A(13)=    8.22733944423959e-06;
	A(14)=      0.999986101965304;
	A(15)=   3.99927483270551e-08;
	A(16)=    3.04673960163327e-06;
	A(17)=   -5.1596062466179e-06;
	A(18)=    1.22075292404534e-08;
	A(19)=   2.17105640392839e-07;
	A(20)=     3.9992870847826e-08;
	A(21)=       0.999984551298692;
	A(22)=    6.63510875632706e-08;
	A(23)=    1.22076480274715e-08;
	A(24)=  -5.73276287738792e-06;
	A(25)=      5.38976081674752;
	A(26)=    1.47507305174403e-05;
	A(27)=   3.21241787851554e-07;
	A(28)=       1.00002936035846;
	A(29)=   8.19365458482084e-06;
	A(30)=    1.40504658112974e-07;
	A(31)=    1.47507306419397e-05;
	A(32)=        5.38968310056198;
	A(33)=    5.90697768748029e-08;
	A(34)=   8.19365482653896e-06;
	A(35)=          0.9999860891763;
	A(36)=    2.58022974647481e-08;
	A(37)=     3.21242427100724e-07;
	A(38)=  5.90698876854246e-08;
	A(39)=         5.38968032557769;
	A(40)=      1.4050537070756e-07;
	A(41)=    2.58024285760964e-08;
	A(42)=        0.999984550703337;
	Matrix R(42);
	R(1) = 573.068082065557;
	R(2) = -3085.15736953138;
	R(3) = -6736.03068347156;
	R(4) = -7.53489822593659;
	R(5) = -1.37294429126638;
	R(6) = -0.0107597986473575;
	R(7) = 1.08925458231315e-05;
	R(8) = 3.04673932160225e-06;
	R(9) = 6.63504292706821e-08;
	R(10) = 2.02239897508587e-06;
	R(11) = 5.61811901849645e-07;
	R(12) = 4.39846387071934e-09;
	R(13) = 3.04673960163327e-06;
	R(14) = -5.1596062466179e-06;
	R(15) = 1.22075292404534e-08;
	R(16) = 5.61812134084449e-07;
	R(17) = -9.58613689243416e-07;
	R(18) = 8.05616500343474e-10;
	R(19) = 6.63510875632706e-08;
	R(20) = 1.22076480274715e-08;
	R(21) = -5.73276287738792e-06;
	R(22) = 4.39895597958216e-09;
	R(23) = 8.0570607835305e-10;
	R(24) = -1.06368693580442e-06;
	R(25) = 1.00002936035846;
	R(26) = 8.19365458482084e-06;
	R(27) = 1.40504658112974e-07;
	R(28) = 1.08999102436198e-05;
	R(29) = 3.02797128053784e-06;
	R(30) = 2.37068516291712e-08;
	R(31) = 8.19365482653896e-06;
	R(32) = 0.9999860891763;
	R(33) = 2.58022974647481e-08;
	R(34) = 3.02797160153579e-06;
	R(35) = -5.16671243316801e-06;
	R(36) = 4.34211426867344e-09;
	R(37) = 1.4050537070756e-07;
	R(38) = 2.58024285760964e-08;
	R(39) = 0.999984550703337;
	R(40) = 2.37075280907946e-08;
	R(41) = 4.34223837651307e-09;
	R(42) = -5.73302112206999e-06;
	R=transpose(R);
	Matrix B = VarEqn( 5.38970808087706,A);

	_assert(m_equals(R,B,abs(R(1)*1e-10)));
    
    return 0;
}
int m_DEInteg_01() {

	Matrix R(6);
	R(1)=5542555.89427451;
	R(2)=3213514.83814162;
	R(3)=3990892.92789074;
	R(4)=5394.06894044389;
	R(5)=-2365.21290574021;
	R(6)=-7061.8448137347;
          
	Matrix A(6);
	A(1)=	6221397.62857869;
	A(2)=	2867713.77965738;
	A(3)=	3006155.98509949;
	A(4)=	4645.04725161806;
	A(5)=  -2752.21591588204;
	A(6)=  -7507.99940987031;
	
	A=transpose(A);
	Matrix B = DEInteg(Accel,0,-134.999991953373,1e-13,1e-6,6,A);


	_assert(m_equals(R,B,abs(R(6)*1e-5)));
    
    return 0;
}
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
    _verify(m_TimeUpdate_01);
    _verify(m_AccelHarmonic_01);
    _verify(m_EqnEquinox_01);
    _verify(m_JPL_Eph_DE430_01);
    _verify(m_LTC_01);
    _verify(m_NutMatrix_01);
    _verify(m_PoleMatrix_01);
    _verify(m_PrecMatrix_01);
    _verify(m_gmst_01);
    _verify(m_gast_01);
    _verify(m_MeasUpdate_01);
    _verify(m_G_AccelHarmonic_01);
    _verify(m_GHAMatrix_01);
    _verify(m_Accel_01); 
    _verify(m_VarEqn_01);
    _verify(m_DEInteg_01);


    return 0;
}


int main()
{
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

	Matrix A(3, 3);
	A(1,1) = 1; A(1,2) = 0; A(1,3) = 0;
	A(2,1) = 2; A(2,2) = 1; A(2,3) = 0;
	A(3,1) = 3; A(3,2) = 14; A(3,3) = 1;

    return (result != 0);
}
