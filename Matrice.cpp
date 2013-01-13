#include "Matrice.hpp"

Matrice::Matrice():m(NULL),N(0){}

Matrice::Matrice(unsigned int N):
	m(new double[N*N]),
	N(N)
{
	std::cout<<"taille"<<std::endl;
	fill_matrice(0.0);
}

Matrice::Matrice(unsigned int N, double val):
	m(new double[N*N]),
	N(N)
{
	std::cout<<"taille+const"<<std::endl;
	fill_matrice(val);
}

Matrice::Matrice(Matrice const& mat):
	m(new double[mat.size()*mat.size()]),
	N(mat.size())
{
	std::cout<<"copie"<<std::endl;
	for(unsigned int i(0);i<N;i++){
		for(unsigned int j(0);j<N;j++){
			m[i*N+j] = mat(i,j);
		}
	}
}

Matrice::~Matrice(){
	delete  m;
	std::cout<<"destructeur"<<std::endl;
}

Matrice& Matrice::operator=(Matrice const& mat){
	std::cout<<"affectation"<<std::endl;
	if(this->N!=mat.N){
		std::cerr<<"impossible d'affecter deux matrices si dim1 =! dim2"<<std::endl;
	}
	for(unsigned int i(0); i<N*N; i++){
		this->m[i] = mat.m[i];
	}
	return (*this);
}

double& Matrice::operator()(unsigned int i, unsigned int j){
	if(i >= N || j>=N){std::cout<<"bug taille"<<std::endl;}
	return m[i*N+j];
}

double const& Matrice::operator()(unsigned int i, unsigned int j) const{
	if(i >= N || j>=N){std::cout<<"bug taille"<<std::endl;}
	return m[i*N+j];
}

Matrice& Matrice::operator*=(Matrice const& mat){
	Matrice tmp(*this);

	for(unsigned int i(0);i<N;i++){
		for(unsigned int j(0);j<N;j++){
			m[i*N+j] = 0.0;
			for(unsigned int k(0);k<N;k++){
				m[i*N+j] += tmp(i,k) * mat(k,j);
			}
		}
	}
	return (*this);
}

Matrice operator*(Matrice const& mat1, Matrice const& mat2){
	Matrice matout(mat1);
	unsigned int N(mat1.size());
	for(unsigned int i(0);i<N;i++){
		for(unsigned int j(0);j<N;j++){
			matout(i,j) = 0.0;
			for(unsigned int k(0);k<N;k++){
				matout(i,j) += mat1(i,k) * mat2(k,j);
			}
		}
	}
	return matout;
}

void Matrice::fill_matrice(double val){
	for(unsigned int i(0);i<N;i++){
		for(unsigned int j(0);j<N;j++){
			m[i*N+j]=val;
		}
	}
}

void Matrice::Test() const{
	for(unsigned int i(0);i<N;i++){
		for(unsigned int j(0);j<N;j++){
			std::cout<<m[i*N+j]<<" ";
		}
		std::cout<<std::endl;;
	}
}

void Matrice::Chop(){
	for(unsigned int i(0);i<N*N;i++){
		if(std::fabs(m[i]) < 1e-10 ){m[i]=0;}
	}
}

Matrice Matrice::Transpose() const{
	Matrice T(N,0.0);
	for(unsigned int i(0);i<N;i++){
		for(unsigned int j(0);j<N;j++){
			T(i,j) = m[j*N+i];
		}
	}
	return T;
}

void Matrice:: eigenvalue(Vecteur& EVal){
	char jobz('V'),uplo('U');
	unsigned int lwork(3*N-1);
	double *work;
	int info;
	work = new double[lwork];
	dsyev_(&jobz, &uplo, &N, m, &N, EVal.ptr(), work ,&lwork, &info);
	delete work;
}

void Matrice::print() const{
	for(unsigned int i(0); i < N; i++){
		for(unsigned int j(0); j < N; j++){
			std::cout << m[i*N+j] << " ";
		}
		std::cout << std::endl;
	}
}

