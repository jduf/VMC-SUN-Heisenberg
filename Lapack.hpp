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

template<typename T>
class Lapack{
	public:
/*Constructors and destructor*/
		Lapack(Matrice<T> const& mat, char matrix_type);
		Lapack(T* m, unsigned int N, char matrix_type);
		~Lapack();

		void eigensystem(Vecteur<double>& EVal, bool EVec=true) const; //by default compute eigenvectors

		T det() const;
		void lu(Matrice<T>& L, Matrice<T>& U) const;
		void inv() const;

		void reset(); 
		
	private:
		Lapack();
		Lapack(Lapack const& l);

		T *m;
		unsigned int const N;
		char const matrix_type; //implémenter les routines lapack selon G,S,H...
		bool const modify_original_matrix;
		
/*private methods used to call lapack*/
		void getrf(int ipiv[]) const;
		void getri(int ipiv[]) const;
		void syev(Vecteur<double> & EVal) const;
		void heev(Vecteur<double> & EVal) const; //valeurs propres d'une matrice hermitienne sont réelles
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
{ }

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
	}
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

template<typename T>
void Lapack<T>::inv() const {
	if (matrix_type != 'G'){
		std::cerr<<"inv : Matrix type "<<matrix_type<<" not implemented"<<std::endl;
		std::cerr<<"inv : the only matrix type implemented is G"<<std::endl;
	}
	int ipiv[N];
	getrf(ipiv);
	getri(ipiv);	
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
