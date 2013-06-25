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
		Matrix(unsigned int total_size);
		/*!Initializes a static array of M of size N*N to a value val*/
		Matrix(unsigned int total_size, Type val);
		/*!Deep copy*/
		Matrix(Matrix const& mat);
		/*!Delete the static array*/
		~Matrix();

		/*!Does a deep copie*/
		//Matrix<Type>& operator=(Matrix<Type> const& mat);
		/*!Accesses the (i,j)th entry of the vector*/
		inline Type const& operator()(unsigned int const& i, unsigned int const& j) const { assert(i<row() && j<col()); return m[i+j*row()]; };
		/*!Sets the (i,j)th entry of the vector*/
		inline Type& operator()(unsigned int const& i, unsigned int const& j) { assert(i<row() && j<col()); return m[i+j*row()]; };

		/*!Set the whole matrix to val*/
		void set(Type const& val);

		/*!Sets the entries to zero if they are close to 0*/
		void chop(double precision = 1e-10);

		/*!Returns the pointer to the matrix*/
		inline Type* ptr() const { return m; };
		/*!Returns the size of the matrix*/
		inline unsigned int size() const { return total_size; };

		/*!Returns the number of rows of the matrix*/
		virtual unsigned int row() const = 0;
		/*!Returns the number of columns of the matrix*/
		virtual unsigned int col() const = 0;

		void test(){std::cout<<"matrix "<<std::endl;}

	protected:
		Type *m; //!< pointer to a static array
		unsigned int total_size; //!< size of the array
};

template<typename Type>
std::ostream& operator<<(std::ostream& flux, Matrix<Type> const& mat);

/*constructors and destructor*/
/*{*/
template<typename Type>
Matrix<Type>::Matrix():
	m(NULL),
	total_size(0)
{ }

template<typename Type>
Matrix<Type>::Matrix(unsigned int total_size):
	m(new Type[total_size]),
	total_size(total_size)
{ } 

template<typename Type>
Matrix<Type>::Matrix(unsigned int total_size, Type val):
	m(new Type[total_size]),
	total_size(total_size)
{
	set(val);
}

template<typename Type>
Matrix<Type>::Matrix(Matrix<Type> const& mat):
	m(new Type[mat.total_size]),
	total_size(total_size)
{
	for(unsigned int i(0);i<total_size;i++){
		m[i] = mat.m[i];
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
std::ostream& operator<<(std::ostream& flux, Matrix<Type> const& mat){
	for(unsigned int i(0);i<mat.row();i++){
		for(unsigned int j(0);j<mat.col();j++){
			flux<<mat(i,j)<<" ";
		}
		flux<<std::endl;
	}
	return flux;
}

//template<typename Type>
//Matrix<Type>& Matrix<Type>::operator=(Matrix<Type> const& mat){
	//if(this->N!=mat.N){
		//if(this->m){ delete[] this->m;}
		//this->m = new Type[mat.N*mat.N];
		//this->N = mat.N;
	//}
	//for(unsigned int i(0); i<this->N*this->N; i++){
		//this->m[i] = mat.m[i];
	//}
	//return (*this);
//}
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


#endif
