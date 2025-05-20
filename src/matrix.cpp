#include "../include/matrix.h"

    /**
     * @file matrix.cpp
     * @brief El archivo contiene las implementaciones de matrix.h
     * @author Pedro Zhuzhan
     * @bug No known bugs
     */
//----------------------------------
Matrix::Matrix() {
	this->n_row = 0;
	this->n_column = 0;
	this->data = NULL;
	
}
//----------------------------------
Matrix::Matrix(const int n_size) {
    if (n_size <= 0) {
		cout << "Vector create: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	this->n_row = 1;
	this->n_column = n_size;
	this->data = (double **) malloc(n_row*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Vector create: error in data\n";
        exit(EXIT_FAILURE);
	}
	this->data[0] = (double *) calloc(n_size,sizeof(double));
	
}
//--------------------------------------
Matrix::Matrix(const int n_row, const int n_column) {
    if (n_row <= 0 || n_column <= 0) {
		cout << "Matrix create: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	this->n_row = n_row;
	this->n_column = n_column;
	this->data = (double **) malloc(n_row*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < n_row; i++) {
		this->data[i] = (double *) malloc(n_column*sizeof(double));
	}
}
//----------------------------------
double& Matrix::operator () (const int n) {
	if (n <= 0 || n > this->n_column* this->n_row) {
		cout << "Vector get: error in get:"<<n<<" size: " <<this->n_column*this->n_row<<" row/column\n";
        exit(EXIT_FAILURE);
	}
	
	return this->data[(n - 1)/this->n_column][(n-1)%this->n_column];
}
//----------------------------------
double& Matrix::operator () (const int row, const int column) {
	if (row <= 0 || row > this->n_row || column <= 0 || column > this->n_column) {
		cout << "Matrix get: error in " <<row<<" row/ "<<column<<" column\n";
		cout << "Matrix " <<this->n_row<<" n_row/ "<<this->n_column<<" n_column\n";
        exit(EXIT_FAILURE);
	}
	
	return this->data[row - 1][column - 1];
}
//----------------------------------
Matrix& Matrix::operator + (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sum: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			m_aux->data[i-1][j-1] = this->data[i-1][j-1] + m.data[i-1][j-1];
		}
	}
	
	return *m_aux;
}
//----------------------------------
Matrix& Matrix::operator - (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sub: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			m_aux->data[i-1][j-1] = this->data[i-1][j-1] - m.data[i-1][j-1];
		}
	}
	
	return *m_aux;
}
//----------------------------------
Matrix& Matrix::operator * (Matrix &m){
	if (this->n_column != m.n_row) {
		cout << "Matrix muliplication: error in n->n_column, m.n_row\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix &m_aux=zeros(this->n_row, m.n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= m.n_column; j++) {
					m_aux.data[i-1][j-1]=0;
        	for(int k = 1; k <= this->n_column; k++) {
				m_aux.data[i-1][j-1] += this->data[i-1][k-1] * m.data[k-1][j-1];
			}
		}
	}
	return m_aux;
}
//----------------------------------
Matrix& Matrix::operator / (Matrix &m){
	if (this->n_column != m.n_row) {
		cout << "Matrix sub: error in n->n_column, m.n_row\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux= new Matrix(this->n_row,m.n_column);
	*m_aux=(*this)*inv(m);
	return *m_aux;
}
//----------------------------------
Matrix& Matrix::operator = (const Matrix &m){
	
	Matrix *m_aux = new Matrix(m.n_row, m.n_column);

	for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++){
					m_aux->data[i-1][j-1]=m.data[i-1][j-1];
        }
    }

	this->n_row = m.n_row;
	this->n_column = m.n_column;
	this->data = (double **) malloc(m.n_row*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < m.n_row; i++) {
		this->data[i] = (double *) malloc(m.n_column*sizeof(double));
	}
	for (int i = 0; i < this->n_row; i++) {
        for (int j = 0; j < this->n_column; j++){
					this->data[i][j]=m_aux->data[i][j];
        }
    }
    return *this;
}
//----------------------------------
Matrix& Matrix::operator + (double d){
    
	Matrix &m_aux=zeros(this->n_row, this->n_column);
	
    for(int i = 0; i < this->n_row; i++) {
        for(int j = 0; j < this->n_column; j++) {
			m_aux.data[i][j] =this->data[i][j]+ d;
		}
	}
	
	return m_aux;
}
Matrix& Matrix::operator - (double d){

	Matrix &m_aux=zeros(this->n_row, this->n_column);

    for(int i = 0; i < this->n_row; i++) {
        for(int j = 0; j < this->n_column; j++) {
			m_aux.data[i][j] =this->data[i][j]- d;
		}
	}
	return m_aux;
}
//----------------------------------
Matrix& Matrix::operator * (double d){
	Matrix &m_aux=zeros(this->n_row, this->n_column);
    for(int i = 0; i < this->n_row; i++) {
        for(int j = 0; j < this->n_column; j++) {
			m_aux.data[i][j] =this->data[i][j]* d;
		}
	}
	return m_aux;
}
//---------------------------------
Matrix& Matrix::operator / (double d){
	Matrix &m_aux=zeros(this->n_row, this->n_column);
    for(int i = 0; i < this->n_row; i++) {
        for(int j = 0; j < this->n_column; j++) {
			m_aux.data[i][j] =this->data[i][j]/ d;
		}
	}
	return m_aux;
}
//----------------------------------
ostream& operator << (ostream &o, Matrix &m) {
	for (int i = 0; i < m.n_row; i++) {
        for (int j = 0; j < m.n_column; j++)
			printf("%5.20lf ", m.data[i][j]);
        o << "\n";
    }
	
    return o;
}
//----------------------------------
Matrix& zeros(const int n_row, const int n_column) {
	Matrix *m_aux = new Matrix(n_row, n_column);
	
	for(int i = 0; i < n_row; i++) {
		for(int j = 0; j < n_column; j++) {
			m_aux->data[i][j] = 0;
		}
	}
	return (*m_aux);
}
//----------------------------------
Matrix& eye(const int size){
	Matrix *m_aux = new Matrix(size,size);
	
	for(int i = 1; i <= size; i++) {
		for(int j = 1; j <= size; j++) {
			m_aux->data[i-1][j-1] = 0;
		}
		m_aux->data[i-1][i-1] = 1;
	}
	
	return (*m_aux);
}
//----------------------------------
Matrix& transpose(Matrix &m) {
	Matrix *m_aux = new Matrix(m.n_column,m.n_row);
	
	for(int i = 0; i < m.n_row; i++) {
		for(int j = 0; j < m.n_column; j++) {
			m_aux->data[j][i] = m.data[i][j];
		}
	}
	
	return (*m_aux);
}
//----------------------------------
void swap_row(Matrix &m,int i,int index){
	double aux;
	for(int k=1;k<=m.n_column;k++){
		aux=m.data[i-1][k-1];
		m.data[i-1][k-1]=m.data[index-1][k-1];
		m.data[index-1][k-1]=aux;
	}
}
//----------------------------------
Matrix& inv(Matrix &m) {
	if (m.n_column != m.n_row) {
		cout << "Matrix sub: error in m.n_column, m.n_row\n";
        exit(EXIT_FAILURE);
	}
	Matrix *m_aux=new Matrix(m.n_row,m.n_column);
	Matrix *m_aux1=new Matrix(m.n_row,m.n_column);
	double ratio,aux;
	*m_aux=m;
	*m_aux1=eye(m.n_row);
	int index;
	for(int i=1;i<=m.n_column;i++){
		index=i;
		aux=fabs(m_aux->data[i-1][i-1]);
		for(int j=i+1;j<=m.n_column;j++){
			if(aux<fabs(m_aux->data[j-1][i-1])){
				aux=fabs(m_aux->data[j-1][i-1]);
				index=j;
			}
		}
		swap_row(*m_aux,i,index);
		swap_row(*m_aux1,i,index);
		if(m_aux->data[i-1][i-1]==0){
			cout << "Error singular Matrix\n";
	        exit(EXIT_FAILURE);
		}
		for(int j=i+1;j<=m.n_column;j++){
			if(m_aux->data[j-1][i-1]!=0){
				ratio=m_aux->data[j-1][i-1]/m_aux->data[i-1][i-1];
				if(ratio!=0){
					for(int k=1;k<m.n_column+1;k++){
						m_aux->data[j-1][k-1]-=ratio*m_aux->data[i-1][k-1];
						m_aux1->data[j-1][k-1]-=ratio*m_aux1->data[i-1][k-1];
					}
				}
			}
		}
	}
	for(int i=m.n_row;i>=1;i--){
		ratio=1/m_aux->data[i-1][i-1];
			for(int k=1;k<=m.n_column;k++){
						m_aux->data[i-1][k-1]*=ratio;
						m_aux1->data[i-1][k-1]*=ratio;
			}
			for(int j=i-1;j>=1;j--){
				ratio=m_aux->data[j-1][i-1];
				for(int k=1;k<=m.n_column;k++){
						m_aux1->data[j-1][k-1]-=ratio*m_aux1->data[i-1][k-1];
				}
			}
		}
		free(m_aux);
	return *m_aux1;
}
//----------------------------------
Matrix& zeros(const int n) {
	Matrix *m_aux = new Matrix(n);
	
	for(int i = 0; i < n; i++) {
			m_aux->data[0][i] = 0;
	}
	
	return (*m_aux);
}
//----------------------------------
double norm(Matrix &m) {
	double r=0;
	
	for(int i = 1; i <= m.n_column*m.n_row; i++) {
			r+=m(i)*m(i);
	}
	return sqrt(r);
}
//----------------------------------
double dot(Matrix &v,Matrix &w){
	if(v.n_column!=w.n_column){
		cout << "Vector dot: error in v.n_column, w.n_column\n";
		exit(EXIT_FAILURE);}
	double result=0;
	for(int i=1;i<=v.n_column*v.n_row;i++)
		result += v(i)*w(i);
	return result;
}

//----------------------------------
Matrix& cross(Matrix &v,Matrix &w){
	if(v.n_column!=w.n_column){
		cout << "Vector cross: error in v.n_column, w.n_column\n";
		exit(EXIT_FAILURE);}
	Matrix *m_aux = new Matrix(v.n_column);
	m_aux->data[0][1] = v(2)*w(3)-w(2)*v(3);
	m_aux->data[0][2] = v(3)*w(1)-w(3)*v(1);
	m_aux->data[0][3] = v(1)*w(2)-w(1)*v(2);
	return (*m_aux);
}
//----------------------------------
    Matrix& extract_vector(Matrix &v,const int start, const int end){
		
	Matrix *m_aux = new Matrix(end-start+1);
	int x=0;
	for (int i=start; i<=end;i++){
		m_aux->data[0][x]=v(i);
		x++;
		}
	return *m_aux;
    }
//----------------------------------
    Matrix& union_vector(Matrix &v,Matrix &w){
    	int x=0,length=v.n_column*v.n_row+w.n_column*w.n_row;
    	Matrix *v_aux=new Matrix(length);
    	for(int i=1; i<=v.n_column*v.n_row;i++){
    		v_aux->data[0][x]=v(i);
    		x++;
    	}
    	for(int i=1; i<=w.n_column*w.n_row;i++){
    		v_aux->data[0][x]=w(i);
    		x++;
    	}
    	return (*v_aux);
    }
//----------------------------------
    Matrix& extract_row(const Matrix &v,const int j){
		if(v.n_row<j || 1>j){
			cout << "Matrix extract_row: error in"<< v.n_row <<" "<<j<<"\n";
			exit(EXIT_FAILURE);}
			Matrix *m_aux = new Matrix(v.n_column);
			for (int i=0;i<v.n_column;i++){
				m_aux->data[0][i]=v.data[j-1][i];
			}
			return (*m_aux);
    }
//----------------------------------
    Matrix& extract_column(Matrix &v,int j){
		if(v.n_column<j || j<1){
			cout << "Matrix extract_column: error in "<< j <<" "<<v.n_column<<"\n";
			exit(EXIT_FAILURE);}
			j--;
			Matrix *m_aux = new Matrix(v.n_row);
			for (int i=0;i<v.n_row;i++){
				m_aux->data[0][i]=v.data[i][j];
			}
			return (*m_aux);
    }
//----------------------------------
    Matrix& assign_row(Matrix &v,Matrix &w,int j){
		if(v.n_row<j || j<1 || v.n_row!=w.n_column){
			cout << "Matrix assign_row: error in v.n_row<j\n";
			exit(EXIT_FAILURE);}
			Matrix &m_aux= v;
			j--;
			for (int i=0;i<m_aux.n_column;i++){
				m_aux.data[j][i]=w(i+1);
			}
			return m_aux;
    }
//----------------------------------
    Matrix& assign_column(Matrix &v,Matrix &w,int j){
		if(v.n_column<j || j<1 || v.n_row!=w.n_column){
			cout << "Matrix assign_column: error in v.n_column<j\n";
			exit(EXIT_FAILURE);}
			Matrix &m_aux= v;
			j--;
			for (int i=0;i<m_aux.n_row;i++){
				m_aux.data[i][j]=w(i+1);
			}
			return m_aux;

    }
//----------------------------------