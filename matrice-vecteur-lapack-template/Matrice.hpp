#ifndef DEF_MATRICE
#define DEF_MATRICE

#include <iostream>
#include <cassert>
#include <cmath> //allow abs(double) and abs(complex) 
#include <complex>

#include "Vecteur.hpp"

/*!Class that implement a static array as a matrix
 *
 * - can be saved with Write.hpp 
 * - can be loaded with Read.hpp 
*/
template<typename Type>
class Matrice{
	public:
		/*!Default constructor that initializes *m to NULL and N to 0*/
		Matrice();
		/*!Initializes a static array of M of size N*N*/
		Matrice(unsigned int N);
		/*!Initializes a static array of M of size N*N to a value val*/
		Matrice(unsigned int N, Type val);
		/*!Deep copy*/
		Matrice(Matrice const& mat);
		/*!Delete the static array*/
		~Matrice();

		/*!Deep copy assignment*/
		Matrice<Type>& operator=(Matrice<Type> const& mat); 
		/*!Accesses the (i,j)th entry of the vector*/
		inline Type const& operator()(unsigned int const& i, unsigned int const& j) const { assert(i<N && j<N); return m[i+j*N]; };
		/*!Sets the (i,j)th entry of the vector*/
		inline Type& operator()(unsigned int const& i, unsigned int const& j) { assert(i<N && j<N); return m[i+j*N]; };
		/*!Additions this matrice with another*/
		Matrice<Type>& operator+=(Matrice<Type> const& mat); 
		Matrice<Type> operator+(Matrice<Type> const& mat) const;
		/*!Substracts this matrice from another (m1 -= m2 : m1 = m1-m2)*/
		Matrice<Type>& operator-=(Matrice<Type> const& mat); 
		Matrice<Type> operator-(Matrice<Type> const& mat) const;
		/*!Multiplies two matrices (m1 *= m2 : m1 = m1*m2)*/
		Matrice<Type>& operator*=(Matrice<Type> const& mat); // 
		Matrice<Type> operator*(Matrice<Type> const& mat) const;

		/*!Sets the entries to zero if they are close to 0*/
		void chop(double precision = 1e-10);

		/*!Set the whole matrix to val*/
		void set(Type const& val);

		/*!Returns the transpose of any matrix*/
		Matrice<Type> transpose() const;
		/*!Returns the conjugate transpose of complex matrix (may give an error) */
		Matrice<std::complex<double> > trans_conj() const;
		/*!Returns the diagonal elements in an vector*/
		Vecteur<Type> diag() const;

		/*!Returns the pointer to the matrix*/
		inline Type* ptr() const { return m; };
		/*!Returns the size of the matrix*/
		inline unsigned int size() const { return N; };

		/*!Print the matrice for mathematica*/
		void print_mathematica();

	private:
		Type *m; //!< pointer to a static array of the form m = [[column0],[column1],...]
		unsigned int N; //!< size of the matrix

};

template<typename Type>
std::ostream& operator<<(std::ostream& flux, Matrice<Type> const& mat);
template<typename Type>
std::istream& operator>>(std::istream& flux, Matrice<Type>& mat);

/*constructors and destructor*/
/*{*/
template<typename Type>
Matrice<Type>::Matrice():
	m(NULL),
	N(0)
{ }

template<typename Type>
Matrice<Type>::Matrice(unsigned int N):
	m(new Type[N*N]),
	N(N)
{ } 

template<typename Type>
Matrice<Type>::Matrice(unsigned int N, Type val):
	m(new Type[N*N]),
	N(N)
{
	set(val);
}

template<typename Type>
Matrice<Type>::Matrice(Matrice<Type> const& mat):
	m(new Type[mat.size()*mat.size()]),
	N(mat.size())
{
	for(unsigned int i(0);i<N*N;i++){
		m[i] = mat.m[i];
	}
}

template<typename Type>
Matrice<Type>::~Matrice(){
	if(m){
		delete[]  m;
		m = NULL;
	}
}
/*}*/

/*operators*/
/*{*/
template<typename Type>
Matrice<Type>& Matrice<Type>::operator=(Matrice<Type> const& mat){
	if(this->N!=mat.N){
		if(this->m){ delete[] this->m;}
		this->m = new Type[mat.N*mat.N];
		this->N = mat.N;
	}
	for(unsigned int i(0); i<this->N*this->N; i++){
		this->m[i] = mat.m[i];
	}
	return (*this);
}

template<typename Type>
Matrice<Type>& Matrice<Type>::operator+=(Matrice<Type> const& mat){
	assert(N == mat.N);
	for(unsigned int i(0);i<N*N;i++){
		this->m[i] += mat.m[i];
	}
	return (*this);
}

template<typename Type>
Matrice<Type> Matrice<Type>::operator+(Matrice<Type> const& mat) const{
	Matrice<Type> matout((*this));
	matout += mat;
	return matout;
}

template<typename Type>
Matrice<Type>& Matrice<Type>::operator-=(Matrice<Type> const& mat){
	assert(N == mat.N);
	for(unsigned int i(0);i<N*N;i++){
		this->m[i] -= mat.m[i];
	}
	return (*this);
}

template<typename Type>
Matrice<Type> Matrice<Type>::operator-(Matrice<Type> const& mat) const{
	Matrice<Type> matout((*this));
	matout -= mat;
	return matout;
}

template<typename Type>
Matrice<Type>& Matrice<Type>::operator*=(Matrice<Type> const& mat){
	assert(N == mat.N);
	Matrice<Type> tmp(*this);
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

template<typename Type>
Matrice<Type> Matrice<Type>::operator*(Matrice<Type> const& mat) const{
	assert(N == mat.N);
	Matrice<Type> matout(N);
	std::cerr<<"Matrice::operator* : l'opération * peut-être améliorée"<<std::endl;
	for(unsigned int i(0);i<N;i++){
		for(unsigned int j(0);j<N;j++){
			Type tmp(0.0);
			for(unsigned int k(0);k<N;k++){
				tmp += this->m[i+k*N] * mat.m[k+j*N];
			}
			matout.m[i+j*N] = tmp;
		}
	}
	return matout;
}

template<typename Type>
std::ostream& operator<<(std::ostream& flux, Matrice<Type> const& mat){
	for(unsigned int i(0);i<mat.size();i++){
		for(unsigned int j(0);j<mat.size();j++){
			flux<<mat(i,j)<<" ";
		}
		flux<<std::endl;
	}
	return flux;
}

template<typename Type>
std::istream& operator>>(std::istream& flux, Matrice<Type>& mat){
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
template<>
inline void Matrice<double>::chop(double precision){
	for(unsigned int i(0);i<N*N;i++){
		if(std::abs(m[i]) < precision ){m[i]=0.0;}
	}
}

template<>
inline void Matrice<std::complex<double> >::chop(double precision){
	for(unsigned int i(0);i<N*N;i++){
		if(std::abs(m[i].imag()) < precision ){m[i].imag()=0.0;}
		if(std::abs(m[i].real()) < precision ){m[i].real()=0.0;}
	}
}

template<typename Type>
void Matrice<Type>::set(Type const& val){
	for(unsigned int i(0); i<this->N*this->N; i++){
		m[i] = val;
	}
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename Type>
Matrice<Type> Matrice<Type>::transpose() const{
	Matrice<Type> tmp(N);
	for(unsigned int i(0);i<N;i++){
		for(unsigned int j(0);j<N;j++){
			tmp.m[i+j*N] = m[j+i*N];
		}
	}
	return tmp;
}

template<typename Type>
Matrice<std::complex<double> > Matrice<Type>::trans_conj() const{
	Matrice<std::complex<double> > tmp(N);
	for(unsigned int i(0);i<N;i++){
		for(unsigned int j(0);j<N;j++){
			tmp.m[i+j*N] = conj(m[j+i*N]);
		}
	}
	return tmp;
}

template<typename Type>
Vecteur<Type> Matrice<Type>::diag() const{
	Vecteur<Type> v(N);
	for(unsigned int i(0);i<N;i++){
		v(i) = m[i*(N+1)];
	}
	return v;
}
/*}*/

template<typename Type>
void Matrice<Type>::print_mathematica(){
	std::cout<<"{{";
	for(unsigned int i(0);i<N-1;i++){
		for(unsigned int j(0);j<N-1;j++){
			std::cout<<m[i+j*N]<<",";
		}
		std::cout<<m[i+(N-1)*N]<<"},"<<std::endl;
		std::cout<<"{";
	}
	for(unsigned int j(0);j<N-1;j++){
		std::cout<<m[N-1+j*N]<<",";
	}
	std::cout<<m[N*N-1]<<"}}"<<std::endl;
}

#endif
