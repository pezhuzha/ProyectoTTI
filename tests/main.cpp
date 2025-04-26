
#include <iomanip>
#include <cmath>
#include "../include/GLOBAL.h"
#include "../include/SAT_Const.h"
#include "../include/IERS.h"
#include "../include/matrix.h"
#include <iostream>

using namespace std;

int main() {
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
	cout<<setprecision(20)<<x_pole<<endl;
	cout<<setprecision(20)<<y_pole<<endl;
	cout<<setprecision(20)<<UT1_UTC<<endl;
	cout<<setprecision(20)<<LOD<<endl;
	cout<<setprecision(20)<<dpsi<<endl;
	cout<<setprecision(20)<<deps<<endl;
	cout<<setprecision(20)<<dx_pole<<endl;
	cout<<setprecision(20)<<dy_pole<<endl;
	cout<<setprecision(20)<<TAI_UTC<<endl;
	cout<<endl;
	cout<<setprecision(20)<<fabs(x_pole-R0)<<endl;
	cout<<setprecision(20)<<fabs(y_pole-R1)<<endl;
	cout<<setprecision(20)<<fabs(UT1_UTC-R2)<<endl;
	cout<<setprecision(20)<<fabs(LOD-R3)<<endl;
	cout<<setprecision(20)<<fabs(dpsi-R4)<<endl;
	cout<<setprecision(20)<<fabs(deps-R5)<<endl;
	cout<<setprecision(20)<<fabs(dx_pole-R6)<<endl;
	cout<<setprecision(20)<<fabs(dy_pole-R7)<<endl;
	cout<<setprecision(20)<<fabs(TAI_UTC-R8)<<endl;
	
    
	
    
    return 0;
}