#include "../include/GLOBAL.h"

    /**
     * @file GLOBAL.cpp
     * @brief El archivo contiene las implementaciones de GLOBAL.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
Matrix eopdata;

void eop19620101(int c){
	eopdata=zeros(13,c);
	
		FILE *fid = fopen("../data/eop19620101.txt","r");
		if(fid==NULL){
			cout << "Fail open eop19620101.txt file \n";
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