#ifndef DEF_MLAPACK
#define DEF_MLAPACK

#include "Matrix.hpp"

//work for a symmetric matrix
/*{*/
extern "C" void dsyev_(
		char const& jobz,
		char const& uplo,
		unsigned int const& n,
		double *m,
		unsigned int const& lda,
		double *w,
		double *work,
		int const& lwork,
		int& info
		);
extern "C" void zheev_(
		char const& jobz,
		char const& uplo,
		unsigned int const& n,
		std::complex<double> *m,
		unsigned int const& lda,
		double *w,
		std::complex<double> *work,
		int const& lwork,
		double *rwork,
		int& info
		);

extern "C" void dgetrf_(
		unsigned int const& row,
		unsigned int const& col,
		double *m,
		unsigned int const& lda,
		int *ipiv,
		int& info
		);
extern "C" void zgetrf_(
		unsigned int const& row,
		unsigned int const& col,
		std::complex<double> *m,
		unsigned int const& lda,
		int *ipiv,
		int& info
		);

extern "C" void dgetri_(
		unsigned int const& n,
		double *m,
		unsigned int const& lda,
		int *ipiv,
		double *work,
		int const& lwork,
		int& info
		);
extern "C" void zgetri_(
		unsigned int const& n,
		std::complex<double> *m,
		unsigned int const& lda,
		int *ipiv,
		std::complex<double> *work,
		int const& lwork,
		int& info
		);

extern "C" void dgecon_(
		char const& norm,
		unsigned int const& n,
		double const *m, 
		unsigned int const& lda,
		double const& anorm, 
		double& rcond, 
		double *work,
		int const *iwork,
		int& info
		);
extern "C" double dlange_(
		char const& norm,
		unsigned int const& row,
		unsigned int const& col,
		double const *m,
		unsigned int const& lda,
		double *work
		);

extern "C" void dgeqp3_(
		unsigned int const& row,
		unsigned int const& col,
		double *m,
		unsigned int const& lda,
		int *jptv,
		double *tau,
		double *work,
		int const& lwork,
		int& info
		);
/*}*/

/*!Class that allows an easy use of the LAPACK routines
 * 
 * - compute determinant
 * - compute inverse
 * - compute eigensystem
 * - compute LU decomposition
 * - compute QR decomposition
 * */
template<typename Type>
class Lapack{
	public:
		/*!Constructor that copy the input matrix, if del=false, the input matrix will be modified*/
		Lapack(Matrix<Type> *m, bool del, char matrix_type);
		/*!Destructor (delete if del==true)*/
		~Lapack();

		/*!Specialized routine that computes the eigenvalues and the eigenvectors if EVec==true*/
		void eigensystem(Matrix<double>& EVal, bool EVec=true); 
		/*!Compute the determinant*/
		Type det();
		/*!Compute the LU decomposition*/
		void lu(Matrix<Type>& L, Matrix<Type>& U);
		/*!Compute the QR decomposition*/
		void qr(Matrix<Type>& Q, Matrix<Type>& R);
		/*!Compute the inverse*/
		void inv();
		/*!Compute the condition number*/
		void cond();
		/*!Compute the norm of the matrix*/
		double norm();

		Matrix<Type>* get_mat(){ return mat; }

		
	protected:
		Matrix<Type> *mat; //!< pointer on a Matrix
		bool del; //!< false if the original matrix will be overwritten by the LAPACK routines
		char const matrix_type; //!< type of matrix (symmetric, hermitian, general,...)
		
		/*!Specialized subroutine that calls a LAPACK routine to compute the LU decomposition*/
		void getrf(int *ipiv);
		/*!Specialized subroutine that calls a LAPACK routine to compute the QR decomposition*/
		void geqp3(Matrix<Type>& Q, Matrix<Type>& R);
		/*!Specialized subroutine that calls a LAPACK routine to compute the inverse after the use of getrf*/
		void getri(int *ipiv);
		/*!Specialized subroutine that calls a LAPACK routine to compute the eigensystem of a symmetric real matrix*/
		void syev(Matrix<double> & EVal) ;
		/*!Specialized subroutine that calls a LAPACK routine to compute the eigensystem of an hermitian complex matrix*/
		void heev(Matrix<double> & EVal) ; 
		/*!Specialized subroutine that calls a LAPACK routine to compute the norm of an general real matrix*/
		double lange() ; 
		/*!Specialized subroutine that calls a LAPACK routine to compute the condition number of an general real matrix*/
		double gecon(double anorm) ; 

		/*!Forbids default constructor*/
		Lapack();
		/*!Forbids copy constructor*/
		Lapack(Lapack const& l);
		/*!Forbids assertion operator*/
		Lapack& operator=(Lapack const& l);

};

/*Constructors and destructor*/
/*{*/
template<typename Type>
Lapack<Type>::Lapack(Matrix<Type> *m, bool del, char matrix_type):
	mat(NULL),
	del(del),
	matrix_type(matrix_type)
{
	if(del){ this->mat = new Matrix<Type>(*m); }
	else { this->mat = m; }
}

template<typename Type>
Lapack<Type>::~Lapack(){
	if(del){ delete this->mat; } 
}
/*}*/

/*methods used to call lapack*/
/*{*/
template<typename Type> 
Type Lapack<Type>::det()  {
	if (matrix_type != 'G'){
		std::cerr<<"Lapack : det : Matrix type "<<matrix_type<<" not implemented"<<std::endl;
		std::cerr<<"Lapack : det : the only matrix type implemented is G"<<std::endl;
		return 0;
	} else {
		/*!
		 * \warning 
		 * the determinant needs to be rewritten for rectangles matrices
		 * */
		Type d(1.0);
		unsigned int N(mat->row());
		int* ipiv(new int[N]);
		getrf(ipiv);
		for(unsigned int i(0); i<N; i++){
			if(ipiv[i] != int(i+1)){ //! \note ipiv return value between [1,N]
				d *= - (*mat)(i,i);
			} else {
				d *= (*mat)(i,i);
			}
		}
		delete[] ipiv;
		return d;
	}
}

template<typename Type>
void Lapack<Type>::inv()  {
	if (matrix_type != 'G'){
		std::cerr<<"Lapack : inv : Matrix type "<<matrix_type<<" not implemented"<<std::endl;
		std::cerr<<"Lapack : inv : the only matrix type implemented is G"<<std::endl;
	} else {
		unsigned int N(mat->row());
		int* ipiv(new int[N]);
		getrf(ipiv);
		getri(ipiv);	
		delete[] ipiv;
	}
}

template<typename Type>
void Lapack<Type>::lu(Matrix<Type>& L, Matrix<Type>& U)  {
	if (matrix_type != 'G'){
		std::cerr<<"Lapack : lu : Matrix type "<<matrix_type<<" not implemented"<<std::endl;
		std::cerr<<"Lapack : lu : the only matrix type implemented is G"<<std::endl;
	} else if (mat->row()!=mat->col()){
		std::cerr<<"Lapack : lu : works for square matrix only"<<std::endl;
	} else  {
		unsigned int N(mat->row());
		L=Matrix<Type>(N,N);
		U=Matrix<Type>(N,N);
		int* ipiv(new int[N]);
		getrf(ipiv);
		for(unsigned int i(0); i<N; i++){
			L(i,i)=1.0;
			U(i,i)=(*mat)[i+i*N];
			for(unsigned int j(i+1); j<N; j++){
				U(i,j) = (*mat)[i+j*N];
				L(j,i) = (*mat)[j+i*N];
			}
		}
		delete[] ipiv;
	}
}

template<typename Type>
void Lapack<Type>::qr(Matrix<Type>& Q, Matrix<Type>& R)  {
	if (matrix_type != 'G'){
		std::cerr<<"Lapack : qr : Matrix type "<<matrix_type<<" not implemented"<<std::endl;
		std::cerr<<"Lapack : qr : the only matrix type implemented is G"<<std::endl;
	} else {
		Matrix<Type> Q12(mat->row(),mat->row(),0.0); //! see wikipedia
		Q = Matrix<Type>(mat->row(),mat->col(),0.0);
		R = Matrix<Type>(mat->col(),mat->col(),0.0);
		geqp3(Q12,R);
		for(unsigned int i(0); i<mat->row(); i++){
			for(unsigned int j(0); j<mat->col(); j++){
				Q(i,j) = Q12(i,j);
			}
		}
	}
}

template<typename Type>
void Lapack<Type>::cond()  {
	if (matrix_type != 'G'){
		std::cerr<<"Lapack : cond : Matrix type "<<matrix_type<<" not implemented"<<std::endl;
		std::cerr<<"Lapack : cond : the only matrix type implemented is G"<<std::endl;
		return 0;
	} else {
		return gecon(lange());
	}
}

template<typename Type> 
double Lapack<Type>::norm()  {
	if (matrix_type != 'G'){
		std::cerr<<"Lapack : norm : Matrix type "<<matrix_type<<" not implemented"<<std::endl;
		std::cerr<<"Lapack : norm : the only matrix type implemented is G"<<std::endl;
		return 0;
	} else {
		return lange();
	}
}
/*}*/
#endif
