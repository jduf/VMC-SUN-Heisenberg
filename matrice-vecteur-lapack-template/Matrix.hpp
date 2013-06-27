#ifndef DEF_MATRIX
#define DEF_MATRIX

#include <iostream>
#include <cassert>
#include <cmath> //allow abs(double) and abs(complex) 
#include <complex>

#include "Vecteur.hpp"

/*!Class that implement a static array as a matrix
 *
 * - can be saved with Write.hpp 
 * - can be loaded with Read.hpp 
*/
template<typename Type>
class Matrix{
	public:
		/*!Default constructor that initializes *m to NULL and N to 0*/
		Matrix();
		/*!Initializes a static array of M of size N*N*/
		Matrix(unsigned int N_row, unsigned int N_col);
		/*!Initializes a static array of M of size N*N to a value val*/
		Matrix(unsigned int N_row, unsigned int N_col, Type val);
		/*!Deep copy*/
		Matrix(Matrix<Type> const& mat);
		/*!Replace a pointer with an instance*/
		//Matrix(Matrix<Type> *mat);
		/*!Delete the static array*/
		~Matrix();

		/*!Accesses the (i,j)th entry of the vector*/
		Type const& operator()(unsigned int const& i, unsigned int const& j) const { assert(i<N_row && j<N_col); return m[i+j*N_row]; };
		/*!Sets the (i,j)th entry of the vector*/
		Type& operator()(unsigned int const& i, unsigned int const& j) { assert(i<N_row && j<N_col); return m[i+j*N_row]; };
		/*!Accesses the (i,j)th entry of the vector*/
		Type const& operator[](unsigned int const& i) const { assert(i<N_total); return m[i]; };
		/*!Sets the (i,j)th entry of the vector*/
		Type& operator[](unsigned int const& i) { assert(i<N_total); return m[i]; };

		/*!Deep copy assignment*/
		Matrix<Type>& operator=(Matrix<Type> const& mat); 
		/*!Additions this matrice with another*/
		Matrix<Type>& operator+=(Matrix<Type> const& mat);
		Matrix<Type> operator+(Matrix<Type> const& mat) const;
		/*!Substracts this matrice from another (m1 -= m2 : m1 = m1-m2)*/
		Matrix<Type>& operator-=(Matrix<Type> const& mat);
		Matrix<Type> operator-(Matrix<Type> const& mat) const;
		/*!Multiplies two matrices (m1 *= m2 : m1 = m1*m2)*/
		Matrix<Type>& operator*=(Matrix<Type> const& mat);
		Matrix<Type> operator*(Matrix<Type> const& mat) const;

		/*!Set the whole matrix to val*/
		void set(Type const& val);
		/*!Sets the entries to zero if they are close to 0*/
		void chop(double precision = 1e-10);
		/*!Print the matrice for mathematica*/
		void print_mathematica();

		/*!Returns the transpose of any matrix*/
		Matrix<Type> transpose() const;
		/*!Returns the conjugate transpose of complex matrix (may give an error) */
		Matrix<Type> trans_conj() const;
		/*!Returns the diagonal elements in an vector*/
		Matrix<Type> diag() const;

		/*!Returns the pointer to the matrix*/
		Type* ptr() const { return m; };
		/*!Returns the size of the matrix*/
		unsigned int size() const { return N_total; };
		/*!Returns the number of rows of the matrix*/
		unsigned int row() const { return N_row; };
		/*!Returns the number of columns of the matrix*/
		unsigned int col() const { return N_col; };

	protected:
		Type *m; //!< pointer to a static array
		unsigned int N_row; //!< number of rows
		unsigned int N_col; //!< number of columns
		unsigned int N_total; //!< size of the array
};

template<typename Type>
std::ostream& operator<<(std::ostream& flux, Matrix<Type> const& mat);
template<typename Type>
std::istream& operator>>(std::istream& flux, Matrix<Type>& mat);

/*constructors and destructor*/
/*{*/
template<typename Type>
Matrix<Type>::Matrix():
	m(NULL),
	N_row(0),
	N_col(0),
	N_total(0)
{
	//std::cout<<"called default"<<std::endl;
}

template<typename Type>
Matrix<Type>::Matrix(unsigned int N_row, unsigned int N_col):
	m(new Type[N_row*N_col]),
	N_row(N_row),
	N_col(N_col),
	N_total(N_col*N_row)
{
	//std::cout<<"called N_row,N_col"<<std::endl;
} 

template<typename Type>
Matrix<Type>::Matrix(unsigned int N_row, unsigned int N_col, Type val):
	m(new Type[N_row*N_col]),
	N_row(N_row),
	N_col(N_col),
	N_total(N_col*N_row)
{
	//std::cout<<"called N_row,N_col,val"<<std::endl;
	set(val);
}

template<typename Type>
Matrix<Type>::Matrix(Matrix<Type> const& mat):
	m(new Type[mat.N_total]),
	N_row(mat.N_row),
	N_col(mat.N_col),
	N_total(mat.N_total)
{
	//std::cout<<"called copy"<<std::endl;
	for(unsigned int i(0);i<N_total;i++){
		this->m[i] = mat.m[i];
	}
}

//template<typename Type>
//Matrix<Type>::Matrix(Matrix<Type> *mat):
	//m(mat->ptr()),
	//N_row(mat->row()),
	//N_col(mat->col()),
	//N_total(mat->size())
//{ 
	//std::cout<<"called pointer"<<std::endl;
//}

template<typename Type>
Matrix<Type>::~Matrix(){
	//std::cout<<"called destructeur";
	if(m){
		//std::cout<<" with delete";
		delete[]  m;
		m = NULL;
	}
	//std::cout<<std::endl;
}
/*}*/

/*operators*/
/*{*/
template<typename Type>
Matrix<Type>& Matrix<Type>::operator=(Matrix<Type> const& mat){
	if(this->N_col != mat.N_col ||  this->N_row != mat.N_row){
		if(this->m){ delete[] this->m;}
		this->m = new Type[mat.N_row*mat.N_col];
		this->N_total = mat.N_total;
		this->N_row = mat.N_row;
		this->N_col = mat.N_col;
	}
	for(unsigned int i(0); i<this->N_total; i++){
		this->m[i] = mat.m[i];
	}
	return (*this);
}

template<typename Type>
Matrix<Type>& Matrix<Type>::operator+=(Matrix<Type> const& mat){
	assert(N_row == mat.N_row && N_col == mat.N_col);
	for(unsigned int i(0);i<N_total;i++){
		this->m[i] += mat.m[i];
	}
	return (*this);
}

template<typename Type>
Matrix<Type> Matrix<Type>::operator+(Matrix<Type> const& mat) const{
	Matrix<Type> matout((*this));
	matout += mat;
	return matout;
}

template<typename Type>
Matrix<Type>& Matrix<Type>::operator-=(Matrix<Type> const& mat){
	assert(N_row == mat.N_row && N_col == mat.N_col);
	for(unsigned int i(0);i<N_total;i++){
		this->m[i] -= mat.m[i];
	}
	return (*this);
}

template<typename Type>
Matrix<Type> Matrix<Type>::operator-(Matrix<Type> const& mat) const{
	Matrix<Type> matout((*this));
	matout -= mat;
	return matout;
}

template<typename Type>
Matrix<Type>& Matrix<Type>::operator*=(Matrix<Type> const& mat){
	assert(this->N_col==mat.N_row);
	Matrix<Type> tmp(*this);
	if(this->N_col!=this->N_row || mat.N_col!=mat.N_row ){
		if(this->m){ delete[] this->m;}
		this->m = new Type[this->N_row*mat.N_col];
		this->N_total = this->N_row*mat.N_col;
		this->N_col = mat.N_col;
	}
	for(unsigned int i(0);i<this->N_row;i++){
		for(unsigned int j(0);j<this->N_col;j++){
			this->m[i+j*this->N_row] = 0.0;
			for(unsigned int k(0);k<tmp.N_col;k++){
				//this->m[i+j*N] += tmp(i,k) * mat(k,j);
				this->m[i+j*this->N_row] += tmp.m[i+k*tmp.N_row] * mat.m[k+j*mat.N_row];
			}
		}
	}
	return (*this);
}

template<typename Type>
Matrix<Type> Matrix<Type>::operator*(Matrix<Type> const& mat) const{
	assert(this->N_col==mat.N_row);
	Matrix<Type> matout(this->N_row,mat.N_col);
	for(unsigned int i(0);i<matout.N_row;i++){
		for(unsigned int j(0);j<matout.N_col;j++){
			matout.m[i+j*matout.N_row] = 0.0;
			for(unsigned int k(0);k<mat.N_row;k++){
				//matout.m[i+j*matout.N_row] += (*this)(i,k) * mat(k,j);
				matout.m[i+j*matout.N_row] += this->m[i+k*this->N_row] * mat.m[k+j*mat.N_row];
			}
		}
	}
	return matout;
}

template<typename Type>
std::ostream& operator<<(std::ostream& flux, Matrix<Type> const& mat){
	for(unsigned int i(0);i<mat.row();i++){
		for(unsigned int j(0);j<mat.col();j++){
			flux<<mat(i,j)<<" ";
		}
		flux<<std::endl;
	}
	return flux;
}

template<typename Type>
std::istream& operator>>(std::istream& flux, Matrix<Type>& mat){
	for(unsigned int i(0);i<mat.row();i++){
		for(unsigned int j(0);j<mat.col();j++){
			flux>>mat(i,j);
		}
	}
	return flux;
}
/*}*/

/*methods that modify the class*/
/*{*/
template<>
inline void Matrix<double>::chop(double precision){
	for(unsigned int i(0);i<N_total;i++){
		if(std::abs(m[i]) < precision ){m[i]=0.0;}
	}
}

template<>
inline void Matrix<std::complex<double> >::chop(double precision){
	for(unsigned int i(0);i<N_total;i++){
		if(std::abs(m[i].imag()) < precision ){m[i].imag()=0.0;}
		if(std::abs(m[i].real()) < precision ){m[i].real()=0.0;}
	}
}

template<typename Type>
void Matrix<Type>::set(Type const& val){
	for(unsigned int i(0); i<N_total; i++){
		m[i] = val;
	}
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename Type>
Matrix<Type> Matrix<Type>::transpose() const{
	Matrix<Type> tmp(this->N_col,this->N_row);
	for(unsigned int i(0);i<N_col;i++){
		for(unsigned int j(0);j<N_row;j++){
			tmp.m[i+j*tmp.N_row] = m[j+i*this->N_row];
		}
	}
	return tmp;
}

template<typename Type>
Matrix<Type> Matrix<Type>::trans_conj() const{
	Matrix<Type> tmp(this->N_col,this->N_row);
	for(unsigned int i(0);i<N_col;i++){
		for(unsigned int j(0);j<N_row;j++){
			tmp.m[i+j*tmp.N_row] = conj(m[j+i*this->N_row]);
		}
	}
	return tmp;
}

template<typename Type>
void Matrix<Type>::print_mathematica(){
	std::cout<<"{{";
	for(unsigned int i(0);i<this->N_row-1;i++){
		for(unsigned int j(0);j<this->N_col-1;j++){
			std::cout<<this->m[i+j*this->N_row]<<",";
		}
		std::cout<<this->m[i+(this->N_col-1)*this->N_row]<<"},"<<std::endl<<"{";
	}
	for(unsigned int j(0);j<this->N_col-1;j++){
		std::cout<<this->m[this->N_row-1+j*this->N_row]<<",";
	}
	std::cout<<this->m[this->N_total-1]<<"}}"<<std::endl;
}

template<typename Type>
Matrix<Type> Matrix<Type>::diag() const{
	unsigned int N(0);
	if(this->N_row < this->N_col){N=N_col;}
	else{N=N_row;}
	Matrix<Type> v(N,1);
	for(unsigned int i(0);i<N;i++){
		v[i] = this->m[i*(this->N_row+1)];
	}
	return v;
}
/*}*/
#endif
