#include "../include/GLOBAL.h"

    /**
     * @file GLOBAL.cpp
     * @brief El archivo contiene las implementaciones de GLOBAL.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
Matrix eopdata;
Matrix Cnm;
Matrix Snm;
Matrix PC;

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

void GGM03S(){
		Cnm=zeros(181,181);
		Snm=zeros(181,181);
		FILE *fid = fopen("../data/GGM03S.txt","r");
		if(fid==NULL){
			cout << "Fail open GGM03S.txt file \n";
			perror("Error");
			exit(EXIT_FAILURE);
		}
		double aux;
		for(int i=1;i<=181;i++){
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
void DE430Coeff(){
	PC=zeros(2285,1020);
		FILE *fid = fopen("../data/DE430Coeff.txt","r");
		if(fid==NULL){
			cout << "Fail open DE430Coeff.txt file \n";
			perror("Error");
			exit(EXIT_FAILURE);
		}
		double aux;
		for(int i=1;i<=2285;i++){
			for (int j=1;j<=1020;j++){
				 fscanf(fid,"%lf",
				 &PC(i,j)
				);
			}
		}
		fclose(fid);
}