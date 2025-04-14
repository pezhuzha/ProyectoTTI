#ifndef _Mjday_
#define _Mjday_

using namespace std;

    /**
     * @file Mjday.h
     * @brief El archivo contiene la funcion Mjday
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
    /**
     * @param year        - year       
     * @param mon         - month      
     * @param day         - day    
     * @param hr          - universal time hour   
     * @param min         - universal time min   
     * @param sec         - universal time sec    
     * @return Modified julian date
     */
	double Mjday(int yr, int mon,int day,int hr=0,int min=0,int sec=0);
#endif


