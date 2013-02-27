#ifndef DEF_MLAPACK
#define DEF_MLAPACK

#include "Matrice.hpp"
#include "Vecteur.hpp"

//work for a symmetric matrix
/*{*/
extern "C" void dsyev_(char *jobz, char *uplo, unsigned int const *n, double *m, unsigned int const *lda, double *w, double work[], unsigned int const *lwork, int *info);
extern "C" void zheev_(char *jobz, char *uplo, unsigned int const *n, std::complex<double> *m, unsigned int const *lda, double *w, std::complex<double> work[], unsigned int const *lwork, double rwork[], int *info);

extern "C" void dgetrf_(unsigned int const *row, unsigned int const *col, double *m, unsigned int const *lda, int ipiv[], int *info);
extern "C" void zgetrf_(unsigned int const *row, unsigned int const *col, std::complex<double> *m, unsigned int const *lda, int ipiv[], int *info);

extern "C" void dgetri_(unsigned int const *n, double *m, unsigned int const *lda, int ipiv[], double work[], unsigned int const *lwork, int *info);
extern "C" void zgetri_(unsigned int const *n, std::complex<double> *m, unsigned int const *lda, int ipiv[], std::complex<double> work[], unsigned int const *lwork, int *info);
/*}*/

/*!Class that allows an easy use of the LAPACK routines
 * 
 * - compute determinant
 * - compute inverse
 * - compute eigensystem
 * - compute LU decomposition
 * */
template<typename T>
class Lapack{
	public:
		/*!Constructor that copy the input matrix, no LAPACK routine can affect the input matrix */
		Lapack(Matrice<T> const& mat, char matrix_type);
		/*!Constructor that takes a pointer on a static array, all the LAPACK routines will affect the static array pointed */
		Lapack(T* m, unsigned int N, char matrix_type);
		/*!Destructor (don't delete anything if modify_original_matrix==false)*/
		~Lapack();

		/*!Specialized routine that computes the eigenvalues and the eigenvectors if EVec==true*/
		void eigensystem(Vecteur<double>& EVal, bool EVec=true) const; 
		/*!Compute the determinant*/
		T det() const;
		/*!Compute the LU decomposition*/
		void lu(Matrice<T>& L, Matrice<T>& U) const;
		/*!Compute the inverse*/
		void inv() const;
		
	private:
		/*!Forbids default constructor*/
		Lapack();
		/*!Forbids copy constructor*/
		Lapack(Lapack const& l);
		/*!Forbids assertion operator*/
		Lapack& operator=(Lapack const& l);

		T *m; //!< pointer on a static array, the square matrix (same structure as the one in Matrice.hpp)
		unsigned int const N; //!< size of the square matrix
		char const matrix_type; //!< type of matrix (symmetric, hermitian, general,...)
		bool const modify_original_matrix; //!< true if the original matrix will be modified by the LAPACK routine
		
		/*!Specialized subroutine that calls a LAPACK routine to compute the LU decomposition*/
		void getrf(int ipiv[]) const;
		/*!Specialized subroutine that calls a LAPACK routine to compute the inverse after the use of getrf*/
		void getri(int ipiv[]) const;
		/*!Specialized subroutine that calls a LAPACK routine to compute the eigensystem of a symmetric real matrix*/
		void syev(Vecteur<double> & EVal) const;
		/*!Specialized subroutine that calls a LAPACK routine to compute the eigensystem of an hermitian complex matrix*/
		void heev(Vecteur<double> & EVal) const; 
};

/*Constructors and destructor*/
/*{*/
template<typename T>
Lapack<T>::Lapack(Matrice<T> const& mat, char matrix_type):
	m(new T[mat.size()*mat.size()]),
	N(mat.size()),
	matrix_type(matrix_type),
	modify_original_matrix(false)
{
	std::cout<<"création (copie matrice) : lapack"<<std::endl;
	std::cout<<"peut-être que la matrice est mal déclarée"<<std::endl;
	for(unsigned int i(0); i<N; i++){
		for(unsigned int j(0); j<N; j++){
			m[i+j*N] = mat(i,j);
		}
	}
	for(unsigned int i(0); i<N; i++){
		for(unsigned int j(0); j<N; j++){
			std::cout<<m[i+j*N]<<" ";
		}
		std::cout<<std::endl;
	}
}

template<typename T>
Lapack<T>::Lapack(T* m, unsigned int N, char matrix_type):
	m(m),
	N(N),
	matrix_type(matrix_type),
	modify_original_matrix(true)
{
	std::cout<<"création (modifie matrice) : lapack"<<std::endl;
}

template<typename T>
Lapack<T>::~Lapack(){
	if(!modify_original_matrix){
		std::cout<<"destructeur : lapack"<<std::endl;
		delete[] m;
	} 
}
/*}*/

/*methods used to call lapack*/
/*{*/
template<typename T> 
T Lapack<T>::det() const {
	if (matrix_type != 'G'){
		std::cerr<<"inv : Matrix type "<<matrix_type<<" not implemented"<<std::endl;
		std::cerr<<"inv : the only matrix type implemented is G"<<std::endl;
		return 0;
	} else {
		T det(1.0);
		int ipiv[N];
		getrf(ipiv);
		for(unsigned int i(0); i<N; i++){
			if(ipiv[i] != i+1){ // ipiv return value between [1,N]
				det *= -m[i*N+i];
			} else {
				det *= m[i*N+i];
			}
		}
		return det;
	}
}

template<typename T>
void Lapack<T>::inv() const {
	if (matrix_type != 'G'){
		std::cerr<<"inv : Matrix type "<<matrix_type<<" not implemented"<<std::endl;
		std::cerr<<"inv : the only matrix type implemented is G"<<std::endl;
	} else {
		int ipiv[N];
		getrf(ipiv);
		getri(ipiv);	
	}
}

template<typename T>
void Lapack<T>::lu(Matrice<T>& L, Matrice<T>& U) const {
	if (matrix_type != 'G'){
		std::cerr<<"inv : Matrix type "<<matrix_type<<" not implemented"<<std::endl;
		std::cerr<<"inv : the only matrix type implemented is G"<<std::endl;
	}
	int ipiv[N];
	getrf(ipiv);
	for(unsigned int i(0); i< N; i++){
		L(i,i)=1.0;
		U(i,i)=m[i+i*N];
		for(unsigned int j(i+1); j< N; j++){
			U(i,j) = m[i+j*N];
			L(j,i) = m[j+i*N];
		}
	}
}
/*}*/

#endif
