
#include <iomanip>
#include <cmath>
#include"../include/matrix.h"
#include"../include/G_AccelHarmonic.h"
#include<tuple>

#include <iostream>

using namespace std;

int main() {
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

	cout<<R-R2<<endl;
	cout<<R<<endl;
	cout<<R2<<endl;

    return 0;
}