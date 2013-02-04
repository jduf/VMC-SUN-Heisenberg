#include "Matrice.hpp"

/*Constructors and destructor*/
/*{*/
Matrice::Matrice():
	m(NULL),
	N(0)
{ } 

Matrice::Matrice(unsigned int N):
	m(new double[N*N]),
	N(N)
{ }

Matrice::Matrice(unsigned int N, double val):
	m(new double[N*N]),
	N(N)
{
	fill_matrice(val);
}

Matrice::Matrice(Matrice const& mat):
	m(new double[mat.size()*mat.size()]),
	N(mat.size())
{
	for(unsigned int i(0);i<N*N;i++){
			m[i] = mat.m[i];
	}
}

Matrice::~Matrice(){
	delete[]  m;
}
/*}*/

/*operators*/
/*{*/
Matrice& Matrice::operator=(Matrice const& mat){
	if(this->N!=mat.N){
		delete[] this->m;
		this->m = new double[mat.N*mat.N];
		this->N = mat.N;
	}
	for(unsigned int i(0); i<this->N*this->N; i++){
		this->m[i] = mat.m[i];
	}
	return (*this);
}

double& Matrice::operator()(unsigned int const &i, unsigned int const &j){
	return m[i+j*N];
}

double const& Matrice::operator()(unsigned int const &i, unsigned int const &j) const{
	return m[i+j*N];
}

Matrice& Matrice::operator+=(Matrice const& mat){
	for(unsigned int i(0);i<N;i++){
		for(unsigned int j(0);j<N;j++){
				m[i+j*N] += mat(i,j);
		}
	}
	return (*this);
}

Matrice operator+(Matrice const& mat1, Matrice const& mat2){
	Matrice matout(mat1);
	matout += mat2;

	return matout;
}

Matrice& Matrice::operator-=(Matrice const& mat){
	for(unsigned int i(0);i<N;i++){
		for(unsigned int j(0);j<N;j++){
				m[i+j*N] -= mat(i,j);
		}
	}
	return (*this);
}

Matrice operator-(Matrice const& mat1, Matrice const& mat2){
	Matrice matout(mat1);
	matout -= mat2;

	return matout;
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

Matrice operator^(Vecteur const& vec1, Vecteur const& vec2){
	if(vec1.size() != vec2.size()){std::cout<<"les vecteurs n'ont pas la même taille"<<std::endl;}
	Matrice mat(vec1.size());
	for(unsigned int i(0);i<mat.size();i++){
		for(unsigned int j(0);j<mat.size();j++){
			mat(i,j) = vec1(i)*vec2(j);
		}
	}
	return mat;
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
			T.m[i+j*N] = m[j+i*N];
		}
	}
	return T;
}
/*}*/

