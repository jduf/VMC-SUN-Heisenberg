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
		Matrix(Matrix const& mat);
		/*!Delete the static array*/
		~Matrix();
		void test();

		/*!Deep copy assignment*/
		Matrix<Type>& operator=(Matrix<Type> const& mat); 
		/*!Accesses the (i,j)th entry of the vector*/
		Type const& operator()(unsigned int const& i, unsigned int const& j) const { assert(i<N_row && j<N_col); return m[i+j*N_row]; };
		/*!Sets the (i,j)th entry of the vector*/
		Type& operator()(unsigned int const& i, unsigned int const& j) { assert(i<N_row && j<N_col); return m[i+j*N_row]; };
		/*!Additions this matrice with another*/
		Matrix<Type>& operator+=(Matrix<Type> const& mat); 
		Matrix<Type> operator+(Matrix<Type> const& mat) const;
		/*!Substracts this matrice from another (m1 -= m2 : m1 = m1-m2)*/
		Matrix<Type>& operator-=(Matrix<Type> const& mat); 
		Matrix<Type> operator-(Matrix<Type> const& mat) const;
		/*!Multiplies two matrices (m1 *= m2 : m1 = m1*m2)*/
		Matrix<Type> operator*(Matrix<Type> const& mat) const;

		/*!Sets the entries to zero if they are close to 0*/
		void chop(double precision = 1e-10);

		/*!Set the whole matrix to val*/
		void set(Type const& val);

		/*!Returns the transpose of any matrix*/
		Matrix<Type> transpose() const;
		/*!Returns the conjugate transpose of complex matrix (may give an error) */
		Matrix<Type> trans_conj() const;

		/*!Returns the pointer to the matrix*/
		Type* ptr() const { return m; };
		/*!Returns the size of the matrix*/
		unsigned int size() const { return total_size; };
		/*!Returns the number of rows of the matrix*/
		unsigned int row() const { return N_row; };
		/*!Returns the number of columns of the matrix*/
		unsigned int col() const { return N_col; };

		/*!Print the matrice for mathematica*/
		void print_mathematica();

	protected:
		Type *m; //!< pointer to a static array
		unsigned int N_row; //!< number of rows
		unsigned int N_col; //!< number of columns
		unsigned int total_size; //!< size of the array
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
	total_size(0)
{ }

template<typename Type>
Matrix<Type>::Matrix(unsigned int N_row, unsigned int N_col):
	m(new Type[N_row*N_col]),
	N_row(N_row),
	N_col(N_col),
	total_size(N_col*N_row)
{ } 

template<typename Type>
Matrix<Type>::Matrix(unsigned int N_row, unsigned int N_col, Type val):
	m(new Type[N_row*N_col]),
	N_row(N_row),
	N_col(N_col),
	total_size(N_col*N_row)
{
	set(val);
}

template<typename Type>
Matrix<Type>::Matrix(Matrix<Type> const& mat):
	m(new Type[mat.total_size]),
	N_row(mat.N_row),
	N_col(mat.N_col),
	total_size(mat.total_size)
{
	for(unsigned int i(0);i<total_size;i++){
		this->m[i] = mat.m[i];
	}
}

template<typename Type>
Matrix<Type>::~Matrix(){
	if(m){
		delete[]  m;
		m = NULL;
	}
}
/*}*/

/*operators*/
/*{*/
template<typename Type>
Matrix<Type>& Matrix<Type>::operator=(Matrix<Type> const& mat){
	if(this->N_col != mat.N_col ||  this->N_row != mat.N_row){
		if(this->m){ delete[] this->m;}
		this->m = new Type[mat.N_row*mat.N_col];
		this->total_size = mat.total_size;
		this->N_row = mat.N_row;
		this->N_col = mat.N_col;
	}
	for(unsigned int i(0); i<this->total_size; i++){
		this->m[i] = mat.m[i];
	}
	return (*this);
}

template<typename Type>
Matrix<Type>& Matrix<Type>::operator+=(Matrix<Type> const& mat){
	assert(N_row == mat.N_row && N_col == mat.N_col);
	for(unsigned int i(0);i<total_size;i++){
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
	for(unsigned int i(0);i<total_size;i++){
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
Matrix<Type> Matrix<Type>::operator*(Matrix<Type> const& mat) const{
	assert(this->N_col==mat.N_row);
	Matrix<Type> matout(this->N_row,mat.N_col);
	for(unsigned int i(0);i<matout.N_row;i++){
		for(unsigned int j(0);j<matout.N_col;j++){
			matout.m[i+j*matout.N_row] = 0.0;
			for(unsigned int k(0);k<mat.N_row;k++){
				//matout.m[i+j*matout.N_row] += (*this)(i,k) * mat(k,j);
				matout.m[i+j*matout.N_row] += m[i+k*this->N_row] * mat.m[k+j*this->N_col];
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
/*}*/

/*methods that modify the class*/
/*{*/
template<>
inline void Matrix<double>::chop(double precision){
	for(unsigned int i(0);i<total_size;i++){
		if(std::abs(m[i]) < precision ){m[i]=0.0;}
	}
}

template<>
inline void Matrix<std::complex<double> >::chop(double precision){
	for(unsigned int i(0);i<total_size;i++){
		if(std::abs(m[i].imag()) < precision ){m[i].imag()=0.0;}
		if(std::abs(m[i].real()) < precision ){m[i].real()=0.0;}
	}
}

template<typename Type>
void Matrix<Type>::set(Type const& val){
	for(unsigned int i(0); i<total_size; i++){
		m[i] = val;
	}
}
/*}*/

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
	for(unsigned int i(0);i<N_row-1;i++){
		for(unsigned int j(0);j<N_col-1;j++){
			std::cout<<m[i+j*N_row]<<",";
		}
		std::cout<<m[i+(N_row-1)*N_row]<<"},"<<std::endl;
		std::cout<<"{";
	}
	for(unsigned int j(0);j<N_col-1;j++){
		std::cout<<m[N_row-1+j*N_row]<<",";
	}
	std::cout<<m[total_size-1]<<"}}"<<std::endl;
}

template<typename Type>
void Matrix<Type>::test(){
	std::cout<<N_row<<" "<<N_col<<" "<<total_size<<std::endl;
}
#endif
