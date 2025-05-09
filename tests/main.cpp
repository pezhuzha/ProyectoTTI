
#include <iomanip>
#include <cmath>
#include"../include/JPL_Eph_DE430.h"
#include"../include/Cheb3D.h"
#include"../include/GLOBAL.h"
#include "../include/SAT_Const.h"
#include "../include/matrix.h"
#include<tuple>

#include <iostream>

using namespace std;

int main() {
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
	
 
cout << "\nMercury:\n" << r_Mercury << "\nReference:\n" << R0 << endl;
cout << "\nVenus:\n"   << r_Venus   << "\nReference:\n" << R1 << endl;
cout << "\nEarth:\n"   << r_Earth   << "\nReference:\n" << R2 << endl;
cout << "\nMars:\n"    << r_Mars    << "\nReference:\n" << R3 << endl;
cout << "\nJupiter:\n" << r_Jupiter << "\nReference:\n" << R4 << endl;
cout << "\nSaturn:\n"  << r_Saturn  << "\nReference:\n" << R5 << endl;
cout << "\nUranus:\n"  << r_Uranus  << "\nReference:\n" << R6 << endl;
cout << "\nNeptune:\n" << r_Neptune << "\nReference:\n" << R7 << endl;
cout << "\nPluto:\n"   << r_Pluto   << "\nReference:\n" << R8 << endl;
cout << "\nMoon:\n"    << r_Moon    << "\nReference:\n" << R9 << endl;
cout << "\nSun:\n"     << r_Sun     << "\nReference:\n" << R10 << endl;

    return 0;
}