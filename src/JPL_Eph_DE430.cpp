#include"../include/JPL_Eph_DE430.h"
#include"../include/Cheb3D.h"
#include"../include/GLOBAL.h"
#include<tuple>

/**
*@file JPL_Eph_DE430.cpp
*@brief El archivo contiene las implementaciones de JPL_Eph_DE430.h
*@author Pedro Zhuzhan
*@bug Noknownbugs
*/
tuple<Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&>	JPL_Eph_DE430(double Mjd_TDB){


double JD,t1,dt,j,Mjd0;
Matrix Cx_Earth(12),Cy_Earth(12),Cz_Earth(12),Cx(12),Cy(12),Cz(12),
Cx_Moon(12),Cy_Moon(12),Cz_Moon(12),temp(4),Cx_Sun(12),Cy_Sun(12),Cz_Sun(12),Cx_Venus(12),Cy_Venus(12),Cz_Venus(12),Cx_Mercury(12),Cy_Mercury(12),Cz_Mercury(12),
Cx_Mars(12),Cy_Mars(12),Cz_Mars(12),Cx_Jupiter(12),Cy_Jupiter(12),Cz_Jupiter(12),Cx_Saturn(12),Cy_Saturn(12),Cz_Saturn(12),Cx_Uranus(12),Cy_Uranus(12),Cz_Uranus(12)
,Cx_Neptune(12),Cy_Neptune(12),Cz_Neptune(12),Cx_Pluto(12),Cy_Pluto(12),Cz_Pluto(12),
&r_Mercury=zeros(3),&r_Venus=zeros(3),&r_Earth=zeros(3),&r_Mars=zeros(3),&r_Jupiter=zeros(3),
&r_Saturn=zeros(3),&r_Uranus=zeros(3),&r_Neptune=zeros(3),&r_Pluto=zeros(3),&r_Moon=zeros(3),&r_Sun=zeros(3);
JD=Mjd_TDB+2400000.5;
int i;
for(i=1;i<=PC.n_row;i++){
	if(PC(i,1)<=JD&&JD<=PC(i,2)){
		break;
	}
}
Matrix PCtemp=extract_row(PC,i);
t1=PCtemp(1)-2400000.5;

dt=Mjd_TDB-t1;
for(int k=1;k<=4;k++){
	temp(k)=231+(k-1)*13;
}
	Cx_Earth=extract_vector(PCtemp,temp(1),temp(2)-1);
	Cy_Earth=extract_vector(PCtemp,temp(2),temp(3)-1);
	Cz_Earth=extract_vector(PCtemp,temp(3),temp(4)-1);
	temp=temp+39;
	Cx=extract_vector(PCtemp,temp(1),temp(2)-1);
	Cy=extract_vector(PCtemp,temp(2),temp(3)-1);
	Cz=extract_vector(PCtemp,temp(3),temp(4)-1);
	Cx_Earth=union_vector(Cx_Earth,Cx);
	Cy_Earth=union_vector(Cy_Earth,Cy);
	Cz_Earth=union_vector(Cz_Earth,Cz);
if(0<=dt&&dt<=16){
	j=0;
Mjd0=t1;
}
else if(16<dt&&dt<=32){
	j=1;
Mjd0=t1+16*j;
}
r_Earth=Cheb3D(Mjd_TDB,13,Mjd0,Mjd0+16,extract_vector(Cx_Earth,13*j+1,13*j+13),extract_vector(Cy_Earth,13*j+1,13*j+13),extract_vector(Cz_Earth,13*j+1,13*j+13))*1e3;


for(int k=1;k<=4;k++){
	temp(k)=441+(k-1)*13;
}
Cx_Moon=extract_vector(PCtemp,temp(1),temp(2)-1);
Cy_Moon=extract_vector(PCtemp,temp(2),temp(3)-1);
Cz_Moon=extract_vector(PCtemp,temp(3),temp(4)-1);
for(i=1;i<7;i++){
	temp=temp+39;
Cx=extract_vector(PCtemp,temp(1),temp(2)-1);
Cy=extract_vector(PCtemp,temp(2),temp(3)-1);
Cz=extract_vector(PCtemp,temp(3),temp(4)-1);
Cx_Moon=union_vector(Cx_Moon,Cx);
Cy_Moon=union_vector(Cy_Moon,Cy);
Cz_Moon=union_vector(Cz_Moon,Cz);
}
if(0<=dt&&dt<=4){
j=0;
Mjd0=t1;}
else if(4<dt&&dt<=8){
j=1;
Mjd0=t1+4*j;}
else if(8<dt&&dt<=12){
j=2;
Mjd0=t1+4*j;}
else if(12<dt&&dt<=16){
j=3;
Mjd0=t1+4*j;}
else if(16<dt&&dt<=20){
j=4;
Mjd0=t1+4*j;}
else if(20<dt&&dt<=24){
j=5;
Mjd0=t1+4*j;}
else if(24<dt&&dt<=28){
j=6;
Mjd0=t1+4*j;}
else if(28<dt&&dt<=32){
j=7;
Mjd0=t1+4*j;}
r_Moon=Cheb3D(Mjd_TDB,13,Mjd0,Mjd0+4,extract_vector(Cx_Moon,13*j+1,13*j+13),extract_vector(Cy_Moon,13*j+1,13*j+13),extract_vector(Cz_Moon,13*j+1,13*j+13))*1e3;
for(int k=1;k<=4;k++){
	temp(k)=753+(k-1)*11;
}
Cx_Sun=extract_vector(PCtemp,temp(1),temp(2)-1);
Cy_Sun=extract_vector(PCtemp,temp(2),temp(3)-1);
Cz_Sun=extract_vector(PCtemp,temp(3),temp(4)-1);
temp=temp+33;
Cx=extract_vector(PCtemp,temp(1),temp(2)-1);
Cy=extract_vector(PCtemp,temp(2),temp(3)-1);
Cz=extract_vector(PCtemp,temp(3),temp(4)-1);
Cx_Sun=union_vector(Cx_Sun,Cx);
Cy_Sun=union_vector(Cy_Sun,Cy);
Cz_Sun=union_vector(Cz_Sun,Cz);
if(0<=dt&&dt<=16){
j=0;
Mjd0=t1;}
else if(16<dt&&dt<=32){
j=1;
Mjd0=t1+16*j;}
r_Sun=Cheb3D(Mjd_TDB,11,Mjd0,Mjd0+16,extract_vector(Cx_Sun,11*j+1,11*j+11),extract_vector(Cy_Sun,11*j+1,11*j+11),extract_vector(Cz_Sun,11*j+1,11*j+11))*1e3;
for(int k=1;k<=4;k++){
	temp(k)=3+(k-1)*14;
}
Cx_Mercury=extract_vector(PCtemp,temp(1),temp(2)-1);
Cy_Mercury=extract_vector(PCtemp,temp(2),temp(3)-1);
Cz_Mercury=extract_vector(PCtemp,temp(3),temp(4)-1);
temp=temp+42;
Cx=extract_vector(PCtemp,temp(1),temp(2)-1);
Cy=extract_vector(PCtemp,temp(2),temp(3)-1);
Cz=extract_vector(PCtemp,temp(3),temp(4)-1);
Cx_Mercury=union_vector(Cx_Mercury,Cx);
Cy_Mercury=union_vector(Cy_Mercury,Cy);
Cz_Mercury=union_vector(Cz_Mercury,Cz);
if(0<=dt&&dt<=8){
j=0;
Mjd0=t1;}
else if(8<dt&&dt<=16){
j=1;
Mjd0=t1+8*j;}
else if(16<dt&&dt<=24){
j=2;
Mjd0=t1+8*j;}
else if(24<dt&&dt<=32){
j=3;
Mjd0=t1+8*j;}
r_Mercury=Cheb3D(Mjd_TDB,14,Mjd0,Mjd0+8,extract_vector(Cx_Mercury,14*j+1,14*j+14),extract_vector(Cy_Mercury,14*j+1,14*j+14),extract_vector(Cz_Mercury,14*j+1,14*j+14))*1e3;
for(int k=1;k<=4;k++){
	temp(k)=171+(k-1)*10;
}
Cx_Venus=extract_vector(PCtemp,temp(1),temp(2)-1);
Cy_Venus=extract_vector(PCtemp,temp(2),temp(3)-1);
Cz_Venus=extract_vector(PCtemp,temp(3),temp(4)-1);
temp=temp+30;
Cx=extract_vector(PCtemp,temp(1),temp(2)-1);
Cy=extract_vector(PCtemp,temp(2),temp(3)-1);
Cz=extract_vector(PCtemp,temp(3),temp(4)-1);
Cx_Venus=union_vector(Cx_Venus,Cx);
Cy_Venus=union_vector(Cy_Venus,Cy);
Cz_Venus=union_vector(Cz_Venus,Cz);
if(0<=dt&&dt<=16){
j=0;
Mjd0=t1;}
else if(16<dt&&dt<=32){
j=1;
Mjd0=t1+16*j;}
r_Venus=Cheb3D(Mjd_TDB,10,Mjd0,Mjd0+16,extract_vector(Cx_Venus,10*j+1,10*j+10),extract_vector(Cy_Venus,10*j+1,10*j+10),extract_vector(Cz_Venus,10*j+1,10*j+10))*1e3;
for(int k=1;k<=4;k++){
	temp(k)=309+(k-1)*11;
}
Cx_Mars=extract_vector(PCtemp,temp(1),temp(2)-1);
Cy_Mars=extract_vector(PCtemp,temp(2),temp(3)-1);
Cz_Mars=extract_vector(PCtemp,temp(3),temp(4)-1);
Mjd0=t1;
r_Mars=Cheb3D(Mjd_TDB,11,Mjd0,Mjd0+32,Cx_Mars,Cy_Mars,Cz_Mars)*1e3;

for(int k=1;k<=4;k++){
	temp(k)=342+(k-1)*8;
}
Cx_Jupiter=extract_vector(PCtemp,temp(1),temp(2)-1);
Cy_Jupiter=extract_vector(PCtemp,temp(2),temp(3)-1);
Cz_Jupiter=extract_vector(PCtemp,temp(3),temp(4)-1);
Mjd0=t1;
r_Jupiter=Cheb3D(Mjd_TDB,8,Mjd0,Mjd0+32,Cx_Jupiter,Cy_Jupiter,Cz_Jupiter)*1e3;

for(int k=1;k<=4;k++){
	temp(k)=366+(k-1)*7;
}
Cx_Saturn=extract_vector(PCtemp,temp(1),temp(2)-1);
Cy_Saturn=extract_vector(PCtemp,temp(2),temp(3)-1);
Cz_Saturn=extract_vector(PCtemp,temp(3),temp(4)-1);
Mjd0=t1;
r_Saturn=Cheb3D(Mjd_TDB,7,Mjd0,Mjd0+32,Cx_Saturn,
Cy_Saturn,Cz_Saturn)*1e3;

for(int k=1;k<=4;k++){
	temp(k)=387+(k-1)*6;
}
Cx_Uranus=extract_vector(PCtemp,temp(1),temp(2)-1);
Cy_Uranus=extract_vector(PCtemp,temp(2),temp(3)-1);
Cz_Uranus=extract_vector(PCtemp,temp(3),temp(4)-1);
Mjd0=t1;
r_Uranus=Cheb3D(Mjd_TDB,6,Mjd0,Mjd0+32,Cx_Uranus,
Cy_Uranus,Cz_Uranus)*1e3;

for(int k=1;k<=4;k++){
	temp(k)=405+(k-1)*6;
}
Cx_Neptune=extract_vector(PCtemp,temp(1),temp(2)-1);
Cy_Neptune=extract_vector(PCtemp,temp(2),temp(3)-1);
Cz_Neptune=extract_vector(PCtemp,temp(3),temp(4)-1);
Mjd0=t1;
r_Neptune=Cheb3D(Mjd_TDB,6,Mjd0,Mjd0+32,Cx_Neptune,
Cy_Neptune,Cz_Neptune)*1e3;
for(int k=1;k<=4;k++){
	temp(k)=423+(k-1)*6;
}
Cx_Pluto=extract_vector(PCtemp,temp(1),temp(2)-1);
Cy_Pluto=extract_vector(PCtemp,temp(2),temp(3)-1);
Cz_Pluto=extract_vector(PCtemp,temp(3),temp(4)-1);
Mjd0=t1;
r_Pluto=Cheb3D(Mjd_TDB,6,Mjd0,Mjd0+32,Cx_Pluto,Cy_Pluto,Cz_Pluto)*1e3;
double EMRAT=81.30056907419062;
double EMRAT1=1/(1+EMRAT);
Matrix aux(3);
aux=r_Moon;
 r_Earth=r_Earth-aux*EMRAT1;
 aux=r_Earth;
 aux=aux*(-1);
 r_Mercury=aux+r_Mercury;
 r_Venus=aux+r_Venus;
 r_Mars=aux+r_Mars;
 r_Jupiter=aux+r_Jupiter;
 r_Saturn=aux+r_Saturn;
 r_Uranus=aux+r_Uranus;
 r_Neptune=aux+r_Neptune;
 r_Pluto=aux+r_Pluto;
 r_Sun=aux+r_Sun;
 return tie(r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,r_Sun);
	}