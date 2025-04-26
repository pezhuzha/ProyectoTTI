
#include <iomanip>
#include <cmath>
#include "../include/NutAngles.h"
#include <iostream>

using namespace std;

int main() {
	double R0 = 2.72256565175042e-05;
	double R1 =  3.87947551912632e-05;

	auto [dpsi, deps]= NutAngles(3);
	cout<<setprecision(20)<<dpsi<<endl;
	cout<<setprecision(20)<<R0<<endl;
	cout<<setprecision(20)<<deps<<endl;
	cout<<setprecision(20)<<R1<<endl;
	cout<<fabs(dpsi-R0)<<endl;
	cout<<fabs(deps-R1);
	
    
    return 0;
}