#include "../include/GLOBAL.h"

    /**
     * @file GLOBAL.cpp
     * @brief El archivo contiene las implementaciones de GLOBAL.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
Param AuxParam;
Matrix eopdata;
Matrix Cnm;
Matrix Snm;
Matrix PC;
void AuxParamLoad(){
	AuxParam.Mjd_UTC=49746.1112847221;
	AuxParam.Mjd_TT=49746.1108586111;
	AuxParam.n=20;
	AuxParam.m=20;
	AuxParam.sun=1;
	AuxParam.moon=1;
	AuxParam.planets=1;
}
void eop19620101(int c){
	eopdata=zeros(13,c);
	
		FILE *fid = fopen("../data/eop19620101.txt","r");
		if(fid==NULL){
			cout << "Fail open eop19620101.txt file \n";
			perror("Error");
			exit(EXIT_FAILURE);
		}
		for (int j=1;j<=c;j++){
			 fscanf(fid,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			&(eopdata(1,j)),&(eopdata (2,j)),&(eopdata (3,j)),
			&(eopdata (4,j)),&(eopdata (5,j)),&(eopdata (6,j)),
			&(eopdata (7,j)),&(eopdata (8,j)),&(eopdata (9,j)),
			&(eopdata (10,j)),&(eopdata (11,j)),&(eopdata (12,j)),
			&(eopdata (13,j))
			);
		}
		fclose(fid);
}

void GGM03S(int n){
		Cnm=zeros(n,n);
		Snm=zeros(n,n);
		FILE *fid = fopen("../data/GGM03S.txt","r");
		if(fid==NULL){
			cout << "Fail open GGM03S.txt file \n";
			perror("Error");
			exit(EXIT_FAILURE);
		}
		double aux;
		for(int i=1;i<=n;i++){
			for (int j=1;j<=i;j++){
				 fscanf(fid,"%lf %lf %lf %lf %lf %lf",
				 &aux,&aux,
				 &Cnm(i,j),&Snm(i,j),
				 &aux,&aux
				);
			}
		}
		fclose(fid);
}
void DE430Coeff(int row,int column){
	PC=zeros(row,column);
		FILE *fid = fopen("../data/DE430Coeff.txt","r");
		if(fid==NULL){
			cout << "Fail open DE430Coeff.txt file \n";
			perror("Error");
			exit(EXIT_FAILURE);
		}
		double aux;
		for(int i=1;i<=row;i++){
			for (int j=1;j<=column;j++){
				 fscanf(fid,"%lf",
				 &PC(i,j)
				);
			}
		}
		fclose(fid);
}