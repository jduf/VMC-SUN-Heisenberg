#ifndef DEF_SQUAREMATRIX
#define DEF_SQUAREMATRIX

#include "Matrix.hpp"

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


		/*!Returns the transpose of any matrix*/
		SquareMatrix<Type> transpose() const;
		/*!Returns the conjugate transpose of complex matrix (may give an error) */
		SquareMatrix<std::complex<double> > trans_conj() const;
		/*!Returns the diagonal elements in an vector*/
		Vecteur<Type> diag() const;

		/*!Returns the number of rows of the matrix*/
		virtual unsigned int row() const { return N;}
		/*!Returns the number of columns of the matrix*/
		virtual	unsigned int col() const { return N;}

		void test(){std::cout<<"square "<<std::endl;}

	private:
		unsigned int N; //!< number of columns/rows in the matrix
};

/*constructors and destructor*/
/*{*/
template<typename Type>
SquareMatrix<Type>::SquareMatrix() { }

template<typename Type>
SquareMatrix<Type>::SquareMatrix(unsigned int N):
	Matrix<Type>(N*N),
	N(N)
{ } 

template<typename Type>
SquareMatrix<Type>::SquareMatrix(unsigned int N, Type val):
	Matrix<Type>(N*N,val),
	N(N)
{ }

template<typename Type>
SquareMatrix<Type>::SquareMatrix(SquareMatrix<Type> const& mat):
	Matrix<Type>(mat),
	N(N)
{ }

template<typename Type>
SquareMatrix<Type>::~SquareMatrix(){ }
/*}*/

#endif
