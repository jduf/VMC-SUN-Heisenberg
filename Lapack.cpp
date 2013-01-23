#include "Lapack.hpp"

/*Constructors and destructor*/
/*{*/
Lapack::Lapack(Matrice const& mat, char matrice_type):
	m(new double[mat.size()*mat.size()]),
	m_original(mat.ptr()),
	N(mat.size()),
	matrice_type(matrice_type),
	preserve_matrix(true)
{
	//std::cout<<"création (copie matrice) : lapack"<<std::endl;
	reset();
}

Lapack::Lapack(double* m, unsigned int N, char matrice_type):
	m(m),
	m_original(m),
	N(N),
	matrice_type(matrice_type),
	preserve_matrix(false)
{
	//std::cout<<"création (modifie matrice) : lapack"<<std::endl;
}

Lapack::~Lapack(){
	if(preserve_matrix){
		//std::cout<<"destructeur : lapack"<<std::endl;
		delete[] m;
	} 
}
/*}*/

/*methods used to call lapack*/
/*{*/
Vecteur Lapack::eigensystem() const {
	if (matrice_type != 'S'){
		std::cerr<<"Matrix type not implemented"<<std::endl;
	}
	return dsyev();
}

double Lapack::det() const {
	int ipiv[N];
	double det(1.0);
	dgetrf(ipiv);
	for(unsigned int i(0); i<N; i++){
		if(ipiv[i] != i+1){
			det *= -m[i*N+i];
		} else {
			det *= m[i*N+i];
		}
	}
	return det;
}

void Lapack::lu(Matrice& L, Matrice& U) const {
	int ipiv[N];
	dgetrf(ipiv);
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

/*private methods used to call lapack*/
/*{*/
void Lapack::dgetrf(int ipiv[]) const {
	int info;
	dgetrf_(&N, &N, m, &N, ipiv, &info);
}

Vecteur Lapack::dsyev() const {
	Vecteur EVal(N);
	char jobz('V'),uplo('U');
	unsigned int const lwork(3*N-1);
	double work[lwork];
	int info;
	dsyev_(&jobz, &uplo, &N, m, &N, EVal.ptr(), work ,&lwork, &info);
	return EVal;
}
/*}*/

/*other methods*/
/*{*/
void Lapack::reset(){
	for(unsigned int i(0); i<N; i++){
		for(unsigned int j(0); j<N; j++){
			m[i+j*N] = m_original[i+j*N];
		}
	}
}

Matrice Lapack::LapToMat() const {
	Matrice mat(N);
	for(unsigned int i(0); i<N; i++){
		for(unsigned int j(0); j<N; j++){
			mat(i,j) = m[i+j*N];
		}
	}

	return mat;
}
/*}*/
