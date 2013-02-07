#ifndef DEF_MATRICE
#define DEF_MATRICE

#include <iostream>
#include <iomanip>
#include <complex>

#include "Vecteur.hpp"

template<typename M>
class Matrice{
	public:
		/*Constructors and destructor*/
		Matrice();
		Matrice(unsigned int N);
		Matrice(unsigned int N, M val);
		Matrice(Matrice const& mat);
		~Matrice();

		/*operators*/
		Matrice<M>& operator=(Matrice<M> const& mat);
		inline M const& operator()(unsigned int const& i, unsigned int const& j) const { return m[i+j*N]; };
		inline M& operator()(unsigned int const& i, unsigned int const& j) { return m[i+j*N]; };
		Matrice<M>& operator+=(Matrice<M> const& mat); 
		Matrice<M>& operator-=(Matrice<M> const& mat); 
		Matrice<M>& operator*=(Matrice<M> const& mat); // m1 *= m2 : m1 = m1*m2

		/*methods that modify the class*/
		void chop();

		/*methods that return something related to the class*/
		Matrice<M> transpose() const;
		Matrice<M> trans_conj() const;

		/*other methods*/
		inline M* ptr() const { return m; };
		inline unsigned int size() const { return N; };

		void print() const;

	private:
		M *m; // m = [[ligne0],[ligne1],...]
		unsigned int N;

		/*methods that modify the class*/
		void fill_matrice(M val);
};

//std::ostream& operator<<(std::ostream& flux, Matrice const& mat);
template<typename M>
Matrice<M> operator+(Matrice<M> const& mat1, Matrice<M> const& mat2);
template<typename M>
Matrice<M> operator-(Matrice<M> const& mat1, Matrice<M> const& mat2);
template<typename M>
Matrice<M> operator*(Matrice<M> const& mat1, Matrice<M> const& mat2);
template<typename M>
Matrice<M> operator^(Vecteur<M> const& vec1, Vecteur<M> const& vec2);


/*Constructors and destructor*/
/*{*/
template<typename M>
Matrice<M>::Matrice():
	m(NULL),
	N(0)
{ }

template<typename M>
Matrice<M>::Matrice(unsigned int N):
	m(new M[N*N]),
	N(N)
{ }

template<typename M>
Matrice<M>::Matrice(unsigned int N, M val):
	m(new M[N*N]),
	N(N)
{ fill_matrice(val); }

template<typename M>
Matrice<M>::Matrice(Matrice<M> const& mat):
	m(new M[mat.size()*mat.size()]),
	N(mat.size())
{
	for(unsigned int i(0);i<N*N;i++){
		m[i] = mat.m[i];
	}
}

template<typename M>
Matrice<M>::~Matrice(){
	delete[]  m;
}
/*}*/

/*operators*/
/*{*/
template<typename M>
Matrice<M>& Matrice<M>::operator=(Matrice<M> const& mat){
	if(this->N!=mat.N){
		delete[] this->m;
		this->m = new M[mat.N*mat.N];
		this->N = mat.N;
	}
	for(unsigned int i(0); i<this->N*this->N; i++){
		this->m[i] = mat.m[i];
	}
	return (*this);
}

template<typename M>
Matrice<M>& Matrice<M>::operator+=(Matrice<M> const& mat){
	for(unsigned int i(0);i<N*N;i++){
		this->m[i] += mat.m[i];
	}
	return (*this);
}

template<typename M>
Matrice<M> operator+(Matrice<M> const& mat1, Matrice<M> const& mat2){
	Matrice<M> matout(mat1);
	matout += mat2;
	return matout;
}

template<typename M>
Matrice<M>& Matrice<M>::operator-=(Matrice<M> const& mat){
	for(unsigned int i(0);i<N*N;i++){
		this->m[i] -= mat.m[i];
	}
	return (*this);
}

template<typename M>
Matrice<M> operator-(Matrice<M> const& mat1, Matrice<M> const& mat2){
	Matrice<M> matout(mat1);
	matout -= mat2;
	return matout;
}

template<typename M>
Matrice<M>& Matrice<M>::operator*=(Matrice<M> const& mat){
	Matrice<M> tmp(*this);

	for(unsigned int i(0);i<N;i++){
		for(unsigned int j(0);j<N;j++){
			this->m[i+j*N] = 0.0;
			for(unsigned int k(0);k<N;k++){
				this->m[i+j*N] += tmp.m[i+k*N] * mat.m[k+j*N];
			}
		}
	}
	return (*this);
}

template<typename M>
Matrice<M> operator*(Matrice<M> const& mat1, Matrice<M> const& mat2){
	Matrice<M> matout(mat1);
	std::cerr<<"opération * peut-être améliorée"<<std::endl;
	//matout *= mat2; // pas sur de quoi est mieux...
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

template<typename M>
Matrice<M> operator^(Vecteur<M> const& vec1, Vecteur<M> const& vec2){
	if(vec1.size() != vec2.size()){std::cout<<"les vecteurs n'ont pas la même taille"<<std::endl;}
	Matrice<M> mat(vec1.size());
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
template<typename M>
void Matrice<M>::chop(){
	std::cerr<<"chop : attention utilise fabs... ne marche peut-être pas avec tous les types"<<std::endl;
	for(unsigned int i(0);i<N*N;i++){
		if(std::fabs(m[i]) < 1e-10 ){m[i]=0;}
	}
}

template<typename M>
void Matrice<M>::fill_matrice(M val){
	for(unsigned int i(0);i<N*N;i++){
		m[i]=val;
	}
}
/*}*/

/*other methods*/
/*{*/
template<typename M>
void Matrice<M>::print() const{
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
template<typename M>
Matrice<M> Matrice<M>::transpose() const{
	Matrice<M> tmp(N);
	for(unsigned int i(0);i<N;i++){
		for(unsigned int j(0);j<N;j++){
			tmp.m[i+j*N] = m[j+i*N];
		}
	}
	return tmp;
}
/*}*/
#endif
