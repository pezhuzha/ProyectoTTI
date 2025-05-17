
#include <iomanip>
#include <cmath>
#include<tuple>
#include <iostream>
#include"../include/matrix.h"
#include "../include/Accel.h"
#include "../include/PrecMatrix.h"
#include "../include/NutMatrix.h"
#include "../include/IERS.h"
#include "../include/timediff.h"
#include "../include/PoleMatrix.h"
#include "../include/AccelHarmonic.h"
#include "../include/GHAMatrix.h"
#include "../include/JPL_Eph_DE430.h"
#include "../include/GLOBAL.h"
#include "../include/AccelPointMass.h"
#include "../include/Mjday_TDB.h"
#include "../include/SAT_Const.h"
#include "../include/Position.h"
#include "../include/Mjday.h"
#include "../include/DEInteg.h"
#include "../include/TimeUpdate.h"
#include "../include/AzElPa.h"
#include "../include/R_x.h"
#include "../include/R_y.h"
#include "../include/R_z.h"
#include "../include/gmst.h"
#include "../include/VarEqn.h"
#include "../include/LTC.h"
#include "../include/MeasUpdate.h"
    /**
     * @file EKF_GEOS3.cpp
     * @brief El archivo contiene las implementaciones de EKF_GEOS3.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
        int main(){

    AuxParamLoad();
    eop19620101();
    GGM03S();
    DE430Coeff();
    GEOS3();

    int i=0,j,ii,nobs = 46;
            double sigma_range,sigma_az,sigma_el,lat,lon,alt,Mjd1,Mjd2,Mjd3,Mjd0,Mjd_UTC=obs(9,1),
            n_eqn,theta,t_old,Mjd_TT,Dist,Mjd_UT1,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC,
            x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,Azim, Elev,t;

            Matrix Rs,Y0_apr,P,LT,yPhi,Phi,Y_true=zeros(6),Y0=zeros(6),U,Y_old,dDdY,dDds,r,s,dEdY,
            K, Y, dAds, dEds,dAdY;
            sigma_range = 92.5;          
sigma_az = 0.0224*Rad; 
sigma_el = 0.0139*Rad; 


lat = Rad*21.5748;     
lon = Rad*(-158.2706); 
alt = 300.20;        
Rs = transpose(Position(lon, lat, alt));

Mjd1 = obs(1,1);
Mjd2 = obs(9,1);
Mjd3 = obs(18,1);
Matrix r2(3),v2(3);
r2(1)=6221397.62857869;
r2(2)=2867713.77965738;
r2(3)=3006155.98509949;

v2(1)= 4645.04725161806;
v2(2)=-2752.21591588204;
v2(3)=-7507.99940987031;
//auto [r2,v2] = anglesg(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);



Y0_apr = union_vector(r2,v2);
        
Mjd0 = Mjday(1995,1,29,02,38,0);

AuxParam.Mjd_UTC = Mjd_UTC;
        
Mjd_UTC = obs(9,1);
n_eqn  = 6;

Y = DEInteg(Accel,0,-(obs(9,1)-Mjd0)*86400.0,1e-13,1e-6,6,Y0_apr);
P = zeros(6,6);
for (i=1;i<=3;i++){
    P(i,i)=1e8;}
for (i=4;i<=6;i++){
    P(i,i)=1e3;}

LT = LTC(lon,lat);

yPhi = zeros(42,1);
Phi  = zeros(6,6);


t = 0;

for (i=1;i<=nobs;i++){
    t_old = t;
    Y_old = Y;
    
    Mjd_UTC = obs(i,1);                       
    t       = (Mjd_UTC-Mjd0)*86400.0;         
    
    tie (x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC) = IERS(eopdata,Mjd_UTC,'l');
    tie (UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC) = timediff(UT1_UTC,TAI_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;
    for ( ii=1;ii<=6;ii++){
        yPhi(ii) = Y_old(ii);
        for (j=1;j<=6;j++)  {
            if (ii==j) {
                yPhi(6*j+ii) = 1; 
            }
            else{
                yPhi(6*j+ii) = 0;
            
            }
        }
    }
    yPhi = DEInteg (VarEqn,0,t-t_old,1e-13,1e-6,42,yPhi);
    
    for (j=1;j<=6;j++){
        Phi = assign_column(Phi,extract_vector(yPhi,6*j+1,6*j+6),j);
    
    }
    Y = DEInteg (Accel,0,t-t_old,1e-13,1e-6,6,Y_old);
    theta = gmst(Mjd_UT1);                    
    U = R_z(theta);
    r = extract_vector(Y,1,3);
    r=transpose(r);
    s = LT*(U*r-Rs); 
  
    P = TimeUpdate(P, Phi);
        
    
    tie( Azim, Elev, dAds, dEds)= AzElPa(s);     
    dAdY = union_vector(dAds*LT*U,zeros(1,3));
    
    
     tie( K, Y, P) = MeasUpdate ( Y, obs(i,2), Azim, sigma_az, dAdY, P, 6 );
    
    r = extract_vector(Y,1,3);
    r=transpose(r);
    s = LT*(U*r-Rs);              
    tie( Azim, Elev, dAds, dEds)= AzElPa(s);


    dEdY = union_vector(dEds*LT*U,zeros(1,3));
    
    
     tie( K, Y, P) = MeasUpdate ( Y, obs(i,3), Elev, sigma_el, dEdY, P, 6 );
    
    r = extract_vector(Y,1,3);
    r=transpose(r);
    s = LT*(U*r-Rs);
                 
    Dist = norm(s);
    dDds = transpose(s/Dist);         
    dDdY = union_vector(dDds*LT*U,zeros(1,3));
    
    
     tie(K, Y,P) = MeasUpdate ( Y, obs(i,4), Dist, sigma_range, dDdY, P, 6 );
}

tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC) = IERS(eopdata,obs(46,1),'l');
tie(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC) = timediff(UT1_UTC,TAI_UTC);
Mjd_TT = Mjd_UTC + TT_UTC/86400;
AuxParam.Mjd_UTC = Mjd_UTC;
AuxParam.Mjd_TT = Mjd_TT;

Y0 = DEInteg (Accel,0,-(obs(46,1)-obs(1,1))*86400.0,1e-13,1e-6,6,Y);

 double aux[]= {5753.173e3, 2673.361e3, 3440.304e3, 4.324207e3, -1.924299e3, -5.728216e3};
 Y_true(6);
 for(i=0;i<6;i++){
     Y_true(i+1)=aux[i];
 }

cout<<"\nError of Position Estimation\n";
cout<<Y0(1)-Y_true(1)<<" [m]\n";
cout<<Y0(2)-Y_true(2)<<" [m]\n";
cout<<Y0(3)-Y_true(3)<<" [m]\n";
cout<<"\nError of Velocity Estimation\n";
cout<<Y0(4)-Y_true(4)<<" [m/s]\n";
cout<<Y0(5)-Y_true(5)<<" [m/s]\n";
cout<<Y0(6)-Y_true(6)<<" [m/s]\n";
return 0;
    }