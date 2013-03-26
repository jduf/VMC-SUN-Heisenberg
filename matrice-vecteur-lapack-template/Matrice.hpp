#ifndef DEF_MATRICE
#define DEF_MATRICE

#include <iostream>
#include <complex>
#include <cassert>

#include "Vecteur.hpp"

/*!Class that implement a static array as a matrix
 *
 * - can be saved with Write.hpp 
 * - can be loaded with Read.hpp 
*/
template<typename M>
class Matrice{
	public:
		/*!Default constructor that initializes *m to NULL and N to 0*/
		Matrice();
		/*!Initializes a static array of M of size N*N*/
		Matrice(unsigned int N);
		/*!Initializes a static array of M of size N*N to a value val*/
		Matrice(unsigned int N, M val);
		/*!Deep copy*/
		Matrice(Matrice const& mat);
		/*!Delete the static array*/
		~Matrice();

		/*!Does a deep copie*/
		Matrice<M>& operator=(Matrice<M> const& mat);
		/*!Accesses the (i,j)th entry of the vector*/
		inline M const& operator()(unsigned int const& i, unsigned int const& j) const { assert(i<N && j<N); return m[i+j*N]; };
		/*!Sets the (i,j)th entry of the vector*/
		inline M& operator()(unsigned int const& i, unsigned int const& j) { assert(i<N && j<N); return m[i+j*N]; };
		/*!Additions this matrice with another*/
		Matrice<M>& operator+=(Matrice<M> const& mat); 
		Matrice<M> operator+(Matrice<M> const& mat) const;
		/*!Substracts this matrice from another (m1 -= m2 : m1 = m1-m2)*/
		Matrice<M>& operator-=(Matrice<M> const& mat); 
		Matrice<M> operator-(Matrice<M> const& mat) const;
		/*!Multiplies two matrices (m1 *= m2 : m1 = m1*m2)*/
		Matrice<M>& operator*=(Matrice<M> const& mat); // 
		Matrice<M> operator*(Matrice<M> const& mat) const;

		/*!Sets the entries to zero if they are close to 0*/
		void chop();

		/*!Set the whole matrix to val*/
		void set(M const& val);

		/*!Returns the transpose of any matrix*/
		Matrice<M> transpose() const;
		/*!Returns the conjugate transpose of complex matrix (may give an error) */
		Matrice<std::complex<double> > trans_conj() const;

		/*!Returns the pointer to the matrix*/
		inline M* ptr() const { return m; };
		/*!Returns the size of the matrix*/
		inline unsigned int size() const { return N; };

	private:
		M *m; //!< pointer to a static array of the form m = [[ligne0],[ligne1],...]
		unsigned int N; //!< size of the matrix

};

template<typename M>
std::ostream& operator<<(std::ostream& flux, Matrice<M> const& mat);
template<typename M>
std::istream& operator>>(std::istream& flux, Matrice<M>& mat);

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
{
	set(val);
}

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
	if(m){
		delete[]  m;
		m = NULL;
	}
}
/*}*/

/*operators*/
/*{*/
template<typename M>
Matrice<M>& Matrice<M>::operator=(Matrice<M> const& mat){
	if(this->N!=mat.N){
		if(this->m){ delete[] this->m;}
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
	assert(N == mat.N);
	for(unsigned int i(0);i<N*N;i++){
		this->m[i] += mat.m[i];
	}
	return (*this);
}

template<typename M>
Matrice<M> Matrice<M>::operator+(Matrice<M> const& mat) const{
	Matrice<M> matout((*this));
	matout += mat;
	return matout;
}

template<typename M>
Matrice<M>& Matrice<M>::operator-=(Matrice<M> const& mat){
	assert(N == mat.N);
	for(unsigned int i(0);i<N*N;i++){
		this->m[i] -= mat.m[i];
	}
	return (*this);
}

template<typename M>
Matrice<M> Matrice<M>::operator-(Matrice<M> const& mat) const{
	Matrice<M> matout((*this));
	matout -= mat;
	return matout;
}

template<typename M>
Matrice<M>& Matrice<M>::operator*=(Matrice<M> const& mat){
	assert(N == mat.N);
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
Matrice<M> Matrice<M>::operator*(Matrice<M> const& mat) const{
	assert(N == mat.N);
	Matrice<M> matout(N,0.0);
	std::cerr<<"Matrice : opération * peut-être améliorée"<<std::endl;
	for(unsigned int i(0);i<N;i++){
		for(unsigned int j(0);j<N;j++){
			for(unsigned int k(0);k<N;k++){
				matout.m[i+j*N] += m[i+k*N] * mat.m[k+j*N];
			}
		}
	}
	return matout;
}

template<typename M>
std::ostream& operator<<(std::ostream& flux, Matrice<M> const& mat){
	for(unsigned int i(0);i<mat.size();i++){
		for(unsigned int j(0);j<mat.size();j++){
			flux<<mat(i,j)<<" ";
		}
		flux<<std::endl;
	}
	return flux;
}

template<typename M>
std::istream& operator>>(std::istream& flux, Matrice<M>& mat){
	for(unsigned int i(0);i<mat.size();i++){
		for(unsigned int j(0);j<mat.size();j++){
			flux>>mat(i,j);
		}
	}
	return flux;
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
void Matrice<M>::set(M const& val){
	for(unsigned int i(0); i<this->N*this->N; i++){
		m[i] = val;
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

template<typename M>
Matrice<std::complex<double> > Matrice<M>::trans_conj() const{
	Matrice<std::complex<double> > tmp(N);
	for(unsigned int i(0);i<N;i++){
		for(unsigned int j(0);j<N;j++){
			tmp.m[i+j*N] = conj(m[j+i*N]);
		}
	}
	return tmp;
}
/*}*/
#endif
