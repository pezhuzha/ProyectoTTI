
#include <iomanip>
#include <cmath>
#include"../include/matrix.h"
#include"../include/Accel.h"
#include<tuple>
#include <iostream>

using namespace std;

int main() {
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
	
	cout<<R<<endl;
	cout<<B<<endl;

    return 0;
}