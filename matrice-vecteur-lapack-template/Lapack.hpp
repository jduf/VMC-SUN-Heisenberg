#ifndef DEF_MLAPACK
#define DEF_MLAPACK

#include "Matrice.hpp"
#include "Vecteur.hpp"

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
		unsigned int const& rows,
		unsigned int const& col,
		double const *m,
		unsigned int const& lda,
		double *work
		);
/*}*/

/*!Class that allows an easy use of the LAPACK routines
 * 
 * - compute determinant
 * - compute inverse
 * - compute eigensystem
 * - compute LU decomposition
 * */
template<typename Type>
class Lapack{
	public:
		/*!Constructor that copy the input matrix, no LAPACK routine can affect the input matrix */
		Lapack(Matrice<Type> const& mat, char matrix_type);
		/*!Constructor that takes a pointer on a static array, all the LAPACK routines will affect the static array pointed */
		Lapack(Type *m, unsigned int N, char matrix_type);
		/*!Destructor (delete if allocate_new_memory==true)*/
		~Lapack();

		/*!Specialized routine that computes the eigenvalues and the eigenvectors if EVec==true*/
		void eigensystem(Vecteur<double>& EVal, bool EVec=true) const; 
		/*!Compute the determinant*/
		Type det() const;
		/*!Compute the LU decomposition*/
		void lu(Matrice<Type>& L, Matrice<Type>& U) const;
		/*!Compute the inverse*/
		void inv() const;
		
		/*!Compute the condition number*/
		void cond() const;

		/*!Compute the norm of the matrix*/
		double norm() const;
		
	private:
		/*!Forbids default constructor*/
		Lapack();
		/*!Forbids copy constructor*/
		Lapack(Lapack const& l);
		/*!Forbids assertion operator*/
		Lapack& operator=(Lapack const& l);

		Type *m; //!< pointer on a static array, the square matrix (same structure as the one in Matrice.hpp)
		unsigned int const N; //!< size of the square matrix
		char const matrix_type; //!< type of matrix (symmetric, hermitian, general,...)
		bool const allocate_new_memory; //!< false if the original matrix will be modified by the LAPACK routine
		
		/*!Specialized subroutine that calls a LAPACK routine to compute the LU decomposition*/
		void getrf(int *ipiv) const;
		/*!Specialized subroutine that calls a LAPACK routine to compute the inverse after the use of getrf*/
		void getri(int *ipiv) const;
		/*!Specialized subroutine that calls a LAPACK routine to compute the eigensystem of a symmetric real matrix*/
		void syev(Vecteur<double> & EVal) const;
		/*!Specialized subroutine that calls a LAPACK routine to compute the eigensystem of an hermitian complex matrix*/
		void heev(Vecteur<double> & EVal) const; 
		
		/*!Specialized subroutine that calls a LAPACK routine to compute the norm of an general real matrix*/
		double lange() const; 

		/*!Specialized subroutine that calls a LAPACK routine to compute the condition number of an general real matrix*/
		double gecon(double anorm) const; 
};

/*Constructors and destructor*/
/*{*/
template<typename Type>
Lapack<Type>::Lapack(Matrice<Type> const& mat, char matrix_type):
	m(new Type[mat.size()*mat.size()]),
	N(mat.size()),
	matrix_type(matrix_type),
	allocate_new_memory(true)
{
	//std::cout<<"création (copie matrice) : lapack"<<std::endl;
	//std::cout<<"peut-être que la matrice est mal déclarée"<<std::endl;
	for(unsigned int i(0); i<N; i++){
		for(unsigned int j(0); j<N; j++){
			m[i+j*N] = mat(i,j);
		}
	}
	//for(unsigned int i(0); i<N; i++){
		//for(unsigned int j(0); j<N; j++){
			//std::cout<<m[i+j*N]<<" ";
		//}
		//std::cout<<std::endl;
	//}
}

template<typename Type>
Lapack<Type>::Lapack(Type *m, unsigned int N, char matrix_type):
	m(m),
	N(N),
	matrix_type(matrix_type),
	allocate_new_memory(false)
{ }

template<typename Type>
Lapack<Type>::~Lapack(){
	if(allocate_new_memory){ delete[] m; } 
}
/*}*/

/*methods used to call lapack*/
/*{*/
template<typename T> 
T Lapack<T>::det() const {
	if (matrix_type != 'G'){
		std::cerr<<"Lapack : det : Matrix type "<<matrix_type<<" not implemented"<<std::endl;
		std::cerr<<"Lapack : det : the only matrix type implemented is G"<<std::endl;
		return 0;
	} else {
		T det(1.0);
		int* ipiv(new int[N]);
		getrf(ipiv);
		for(unsigned int i(0); i<N; i++){
			if(ipiv[i] != int(i+1)){ // ipiv return value between [1,N]
				det *= -m[i*N+i];
			} else {
				det *= m[i*N+i];
			}
		}
		delete[] ipiv;
		return det;
	}
}

template<typename T>
void Lapack<T>::inv() const {
	if (matrix_type != 'G'){
		std::cerr<<"Lapack : inv : Matrix type "<<matrix_type<<" not implemented"<<std::endl;
		std::cerr<<"Lapack : inv : the only matrix type implemented is G"<<std::endl;
	} else {
		int* ipiv(new int[N]);
		getrf(ipiv);
		getri(ipiv);	
		delete[] ipiv;
	}
}

template<typename T>
void Lapack<T>::lu(Matrice<T>& L, Matrice<T>& U) const {
	if (matrix_type != 'G'){
		std::cerr<<"Lapack : lu : Matrix type "<<matrix_type<<" not implemented"<<std::endl;
		std::cerr<<"Lapack : lu : the only matrix type implemented is G"<<std::endl;
	} else {
		int* ipiv(new int[N]);
		getrf(ipiv);
		for(unsigned int i(0); i< N; i++){
			L(i,i)=1.0;
			U(i,i)=m[i+i*N];
			for(unsigned int j(i+1); j< N; j++){
				U(i,j) = m[i+j*N];
				L(j,i) = m[j+i*N];
			}
		}
		delete[] ipiv;
	}
}

template<typename T>
void Lapack<T>::cond() const {
	if (matrix_type != 'G'){
		std::cerr<<"Lapack : cond : Matrix type "<<matrix_type<<" not implemented"<<std::endl;
		std::cerr<<"Lapack : cond : the only matrix type implemented is G"<<std::endl;
		return 0;
	} else {
		return gecon(lange());
	}
}

template<typename T> 
double Lapack<T>::norm() const {
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
