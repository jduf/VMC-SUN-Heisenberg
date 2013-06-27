#ifndef DEF_SQUAREMATRIX
#define DEF_SQUAREMATRIX

#include "Matrix.hpp"
#include "SquareMatrix.hpp"

/*!Class that implement a static array as a matrix
 *
 * - can be saved with Write.hpp 
 * - can be loaded with Read.hpp 
*/
template<typename Type>
class SquareMatrix : public Matrix<Type>{
	public:
		/*!Default constructor that initializes *m to NULL and N to 0*/
		SquareMatrix();
		/*!Initializes a static array of M of size N*N*/
		SquareMatrix(unsigned int N);
		/*!Initializes a static array of M of size N*N to a value val*/
		SquareMatrix(unsigned int N, Type val);
		/*!Deep copy*/
		SquareMatrix(SquareMatrix const& mat);
		/*!Delete the static array*/
		~SquareMatrix();

		/*!Multiplies two matrices (m1 *= m2 : m1 = m1*m2)*/
		SquareMatrix<Type>& operator*=(SquareMatrix<Type> const& mat);  


	private:
		unsigned int N; //!< number of columns/rows in the matrix
};

/*constructors and destructor*/
/*{*/
template<typename Type>
SquareMatrix<Type>::SquareMatrix() { }

template<typename Type>
SquareMatrix<Type>::SquareMatrix(unsigned int N):
	Matrix<Type>(N,N),
	N(N)
{ } 

template<typename Type>
SquareMatrix<Type>::SquareMatrix(unsigned int N, Type val):
	Matrix<Type>(N,N,val),
	N(N)
{ }

template<typename Type>
SquareMatrix<Type>::SquareMatrix(SquareMatrix<Type> const& mat):
	Matrix<Type>(mat),
	N(mat.N)
{ }

template<typename Type>
SquareMatrix<Type>::~SquareMatrix(){ }
/*}*/

/*operators*/
/*{*/
template<typename Type>
SquareMatrix<Type>& SquareMatrix<Type>::operator*=(SquareMatrix<Type> const& mat){
	assert(this->N == mat.N);
	SquareMatrix<Type> tmp(*this);
	for(unsigned int i(0);i<N;i++){
		for(unsigned int j(0);j<N;j++){
			this->m[i+j*N] = 0.0;
			for(unsigned int k(0);k<N;k++){
				//this->m[i+j*N] += tmp(i,k) * mat(k,j);
				this->m[i+j*N] += tmp.m[i+k*N] * mat.m[k+j*N];
			}
		}
	}
	return (*this);
}

/*}*/
#endif
