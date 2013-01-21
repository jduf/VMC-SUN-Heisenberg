#include "Matrice.hpp"

/*Constructors and destructor*/
/*{*/
Matrice::Matrice():
	m(NULL),
	N(0)
{
	std::cout<<"default : matrice"<<std::endl;
	//fill_matrice(0.0);
}

Matrice::Matrice(unsigned int N):
	m(new double[N*N]),
	N(N)
{
	std::cout<<"taille : matrice"<<std::endl;
	//fill_matrice(0.0);
}

Matrice::Matrice(unsigned int N, double val):
	m(new double[N*N]),
	N(N)
{
	std::cout<<"taille+const : matrice"<<std::endl;
	fill_matrice(val);
}

Matrice::Matrice(Matrice const& mat):
	m(new double[mat.size()*mat.size()]),
	N(mat.size())
{
	std::cout<<"copie : matrice"<<std::endl;
	for(unsigned int i(0);i<N;i++){
		for(unsigned int j(0);j<N;j++){
			m[i+j*N] = mat(i,j);
		}
	}
}

Matrice::~Matrice(){
	delete  m;
	std::cout<<"destructeur : matrice"<<std::endl;
}
/*}*/

/*operators*/
/*{*/
Matrice& Matrice::operator=(Matrice const& mat){
	std::cout<<"affectation"<<std::endl;
	if(this->N!=mat.N){
		delete this->m;
		this->m = new double[mat.N*mat.N];
		this->N = mat.N;
		std::cerr<<"affectation d'une matrice de taille différente"<<std::endl;
	}
	for(unsigned int i(0); i<this->N*this->N; i++){
		this->m[i] = mat.m[i];
	}
	return (*this);
}

double& Matrice::operator()(unsigned int i, unsigned int j){
	if(i >= N || j>=N){ std::cout<<"bug taille"<<std::endl; }
	return m[i+j*N];
}

double const& Matrice::operator()(unsigned int i, unsigned int j) const{
	if(i >= N || j>=N){std::cout<<"bug taille"<<std::endl;}
	return m[i+j*N];
}

Matrice& Matrice::operator*=(Matrice const& mat){
	Matrice tmp(*this);

	for(unsigned int i(0);i<N;i++){
		for(unsigned int j(0);j<N;j++){
			m[i+j*N] = 0.0;
			for(unsigned int k(0);k<N;k++){
				m[i+j*N] += tmp(i,k) * mat(k,j);
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
/*}*/

/*methods that modify the class*/
/*{*/
void Matrice::chop(){
	for(unsigned int i(0);i<N*N;i++){
		if(std::fabs(m[i]) < 1e-10 ){m[i]=0;}
	}
}

void Matrice::fill_matrice(double val){
	for(unsigned int i(0);i<N*N;i++){
		m[i]=val;
	}
}
/*}*/

/*other methods*/
/*{*/
void Matrice::print() const{
	for(unsigned int i(0); i < N; i++){
		for(unsigned int j(0); j < N; j++){
			std::cout <<std::setprecision(3)<<std::setw(7)<<std::fixed<< m[i+j*N] << " ";
		}
		std::cout << std::endl;
	}
}
/*}*/

/*methods that return something related to the class*/
/*{*/
Matrice Matrice::transpose() const{
	Matrice T(N);
	for(unsigned int i(0);i<N;i++){
		for(unsigned int j(0);j<N;j++){
			T(i,j) = m[j+i*N];
		}
	}
	return T;
}
/*}*/

