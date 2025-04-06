#include "../include/matrix.h"
//----------------------------------
Matrix::Matrix(const int n_size) {
    if (n_size <= 0) {
		cout << "Vector create: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	this->n_row = 1;
	this->n_column = n_column;
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
		cout << "Vector get: error in row/column\n";
        exit(EXIT_FAILURE);
	}
	
	return this->data[(n - 1)/this->n_column][(n-1)%this->n_column];
}
//----------------------------------
double& Matrix::operator () (const int row, const int column) {
	if (row <= 0 || row > this->n_row || column <= 0 || column > this->n_column) {
		cout << "Matrix get: error in row/column\n";
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
			(*m_aux)(i,j) = (*this)(i,j) + m(i,j);
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
			(*m_aux)(i,j) = (*this)(i,j) - m(i,j);
		}
	}
	
	return *m_aux;
}
//----------------------------------
Matrix& Matrix::operator * (Matrix &m){
	if (this->n_column != m.n_row) {
		cout << "Matrix sub: error in n->n_column, m.n_row\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, m.n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= m.n_column; j++) {
					(*m_aux)(i,j)=0;
        	for(int k = 1; k <= this->n_column; k++) {
			(*m_aux)(i,j) += (*this)(i,k) * m(k,j);
			}
		}
	}
	return *m_aux;
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
Matrix& Matrix::operator = (Matrix &m){

	for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++){
					(*this)(i,j)=m(i,j);
        }
    }
    return *this;
}
//----------------------------------
Matrix& Matrix::operator + (double d){
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*this)(i,j) += d;
		}
	}
	return *this;
}
Matrix& Matrix::operator - (double d){
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*this)(i,j) -= d;
		}
	}
	return *this;
}
//----------------------------------
Matrix& Matrix::operator * (double d){
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*this)(i,j) *= d;
		}
	}
	return *this;
}
//---------------------------------
Matrix& Matrix::operator / (double d){
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*this)(i,j) /= d;
		}
	}
	return *this;
}
//----------------------------------
ostream& operator << (ostream &o, Matrix &m) {
	for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++)
			printf("%5.20lf ", m(i,j));
        o << "\n";
    }
	
    return o;
}
//----------------------------------
Matrix& zeros(const int n_row, const int n_column) {
	Matrix *m_aux = new Matrix(n_row, n_column);
	
	for(int i = 1; i <= n_row; i++) {
		for(int j = 1; j <= n_column; j++) {
			(*m_aux)(i,j) = 0;
		}
	}
	
	return (*m_aux);
}
//----------------------------------
Matrix& eye(const int size){
	Matrix *m_aux = new Matrix(size,size);
	
	for(int i = 1; i <= size; i++) {
		for(int j = 1; j <= size; j++) {
			(*m_aux)(i,j) = 0;
		}
		(*m_aux)(i,i) = 1;
	}
	
	return (*m_aux);
}
//----------------------------------
Matrix& transpose(Matrix &m) {
	Matrix *m_aux = new Matrix(m.n_column,m.n_row);
	
	for(int i = 1; i <= m.n_row; i++) {
		for(int j = 1; j <= m.n_column; j++) {
			(*m_aux)(j,i) = m(i,j);
		}
	}
	
	return (*m_aux);
}
//----------------------------------
void swap_row(Matrix &m,int i,int index){
	double aux;
	for(int k=1;k<=m.n_column;k++){
		aux=m(i,k);
		m(i,k)=m(index,k);
		m(index,k)=aux;
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
		aux=fabs((*m_aux)(i,i));
		for(int j=i+1;j<=m.n_column;j++){
			if(aux<fabs((*m_aux)(j,i))){
				aux=fabs((*m_aux)(j,i));
				index=j;
			}
		}
		swap_row(*m_aux,i,index);
		swap_row(*m_aux1,i,index);
		if((*m_aux)(i,i)==0){
			cout << "Error singular Matrix\n";
	        exit(EXIT_FAILURE);
		}
		for(int j=i+1;j<=m.n_column;j++){
			if((*m_aux)(j,i)!=0){
				ratio=(*m_aux)(j,i)/(*m_aux)(i,i);
				if(ratio!=0){
					for(int k=1;k<m.n_column+1;k++){
						(*m_aux)(j,k)-=ratio*(*m_aux)(i,k);
						(*m_aux1)(j,k)-=ratio*(*m_aux1)(i,k);
					}
				}
			}
		}
	}
	for(int i=m.n_row;i>=1;i--){
		ratio=1/(*m_aux)(i,i);
			for(int k=1;k<=m.n_column;k++){
						(*m_aux)(i,k)*=ratio;
						(*m_aux1)(i,k)*=ratio;
			}
			for(int j=i-1;j>=1;j--){
				ratio=(*m_aux)(j,i);
				for(int k=1;k<=m.n_column;k++){
						(*m_aux)(j,k)-=ratio*(*m_aux)(i,k);
						(*m_aux1)(j,k)-=ratio*(*m_aux1)(i,k);
				}
			}
		}
		free(m_aux);
	return *m_aux1;
}
//----------------------------------
Matrix& zeros(const int n) {
	Matrix *m_aux = new Matrix(n);
	
	for(int i = 1; i <= n; i++) {
			(*m_aux)(1,i) = 0;
	}
	
	return (*m_aux);
}
//----------------------------------
double norm(Matrix &m) {
	double r=0;
	
	for(int i = 1; i <= m.n_column; i++) {
			r+=m(1,i)*m(1,i);
	}
	return sqrt(r);
}
//----------------------------------
double dot(Matrix &v,Matrix &w){
	if(v.n_column!=w.n_column)
		exit(EXIT_FAILURE);
	double result=0;
	for(int i=1;i<v.n_column;i++)
		result += v(1,i)*w(1,i);
	return result;
}

//----------------------------------
Matrix& cross(Matrix &v,Matrix &w){
	if(v.n_column!=w.n_column)
		exit(EXIT_FAILURE);
	Matrix *m_aux = new Matrix(v.n_column);
	(*m_aux)(1) = v(2)*w(3)-w(2)*v(3);
	(*m_aux)(2) = v(3)*w(1)-w(3)*v(1);
	(*m_aux)(3) = v(1)*w(2)-w(1)*v(2);
	return (*m_aux);
}