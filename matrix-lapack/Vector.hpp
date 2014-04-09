#ifndef DEF_VECTOR
#define DEF_VECTOR

#include <iostream>
#include <cassert>
#include <cmath> //allow abs(double) and abs(complex) 
#include <complex>

//#include "Matrix.hpp"

template<typename Type>
class Matrix;

/*{!Class that implement a static array as a Vector
 *
 * - can be saved with Write.hpp 
 * - can be loaded with Read.hpp 
 }*/
template<typename Type>
class Vector{
	public:
		/*!Default constructor that initializes *m to NULL and N to 0*/
		Vector();
		/*!Initializes a static array of Type of size N*/
		Vector(unsigned int N);
		/*!Initializes a static array of Type of size N to a value val*/
		Vector(unsigned int N, Type val);
		/*!Deep copy*/
		Vector(Vector<Type> const& vec);
		/*!Delete the static array*/
		~Vector();

		/*!Accesses the (i,j)th entry of the Vector*/
		Type const& operator()(unsigned int const& i) const { 
			assert(i<N_); return m_[i]; };
		/*!Sets the (i,j)th entry of the Vector*/
		Type& operator()(unsigned int const& i) { 
			assert(i<N_); return m_[i]; };

		/*!Deep copy assignment*/
		Vector<Type>& operator=(Vector<Type> const& vec); 
		/*!Additions this vecrice with another*/
		Vector<Type>& operator+=(Vector<Type> const& vec);
		Vector<Type> operator+(Vector<Type> const& vec) const;
		/*!Substracts this vecrice from another (m1 -= m2 : m1 = m1-m2)*/
		Vector<Type>& operator-=(Vector<Type> const& vec);
		Vector<Type> operator-(Vector<Type> const& vec) const;
		/*!Multiplies two vectors (v1 *= v2 : v1 = v1*v2)*/
		Vector<Type> operator*(Vector<Type> const& vec) const;
		/*!Multiplies two vectors and get a matrix*/
		Matrix<Type> operator^(Vector<Type> const& vec) const;

		/*!Devides a vectors by a scalar*/
		Vector<Type>& operator/=(Type const& d);
		Vector<Type> operator/(Type const& d) const;

		/*!Set the whole Vector to val*/
		void set(unsigned int N, Type const& val);
		/*!Set the whole Vector to val*/
		void set(unsigned int N);
		void set();
		/*!Sets the entries to zero if they are close to 0*/
		Vector<Type> chop(double precision = 1e-10) const;
		/*!Print the vector for vechevecica*/
		void print_mathematica();

		Vector<Type> range(unsigned int min, unsigned int max) const;
		Vector<unsigned int> sort();
		Vector<Type> sort(Vector<unsigned int> const& index) const;
		bool is_sorted() const;
		Type mean() const;
		Type max() const;
		Type min() const;

		/*!Returns the pointer to the Vector*/
		Type* ptr() const { return m_; }
		/*!Returns the size of the Vector*/
		unsigned int size() const { return N_; }

	protected:
		Type *m_; //!< pointer to a static array
		unsigned int N_; //!< number of rows

		void set_null_pointer(){m_=NULL;}
};

template<typename Type>
std::ostream& operator<<(std::ostream& flux, Vector<Type> const& vec);
template<typename Type>
std::istream& operator>>(std::istream& flux, Vector<Type>& vec);

/*constructors and destructor*/
/*{*/
template<typename Type>
Vector<Type>::Vector():
	m_(NULL),
	N_(0)
{ }

template<typename Type>
Vector<Type>::Vector(unsigned int N):
	m_(new Type[N]),
	N_(N)
{ } 

template<typename Type>
Vector<Type>::Vector(unsigned int N, Type val):
	m_(new Type[N]),
	N_(N)
{
	set(N,val);
}

template<typename Type>
Vector<Type>::Vector(Vector<Type> const& vec):
	m_(new Type[vec.N_]),
	N_(vec.N_)
{
	for(unsigned int i(0);i<N_;i++){
		m_[i] = vec.m_[i];
	}
}

template<typename Type>
Vector<Type>::~Vector(){
	if(m_){
		delete[]  m_;
		m_ = NULL;
	}
}
/*}*/

/*operators*/
/*{*/
template<typename Type>
Vector<Type>& Vector<Type>::operator=(Vector<Type> const& vec){
	if(!vec.m_){
		set();
	} else {
		if(N_ != vec.N_){
			if(m_){ delete[] m_;}
			m_ = new Type[vec.N_];
			N_ = vec.N_;
		}
		for(unsigned int i(0); i<N_; i++){
			m_[i] = vec.m_[i];
		}
	}
	return (*this);
}

template<typename Type>
Vector<Type>& Vector<Type>::operator+=(Vector<Type> const& vec){
	assert(N_ == vec.N_);
	for(unsigned int i(0);i<N_;i++){
		m_[i] += vec.m_[i];
	}
	return (*this);
}

template<typename Type>
Vector<Type> Vector<Type>::operator+(Vector<Type> const& vec) const{
	Vector<Type> vecout((*this));
	vecout += vec;
	return vecout;
}

template<typename Type>
Vector<Type>& Vector<Type>::operator-=(Vector<Type> const& vec){
	assert(N_ == vec.N_);
	for(unsigned int i(0);i<N_;i++){
		m_[i] -= vec.m_[i];
	}
	return (*this);
}

template<typename Type>
Vector<Type> Vector<Type>::operator-(Vector<Type> const& vec) const{
	Vector<Type> vecout((*this));
	vecout -= vec;
	return vecout;
}

template<typename Type>
Vector<Type> Vector<Type>::operator*(Vector<Type> const& vec) const{
	Type out(0.0);
	for(unsigned int i(0);i<N_;i++){
		out += m_[i] * vec.m_[i];
	}
	return out;
}

template<typename Type>
Matrix<Type> Vector<Type>::operator^(Vector<Type> const& vec) const{
	Matrix<Type> out(vec.size(),vec.size());
	for(unsigned int i(0);i<N_;i++){
		for(unsigned int j(0);j<N_;j++){
			out(i,j) = m_[i] * vec.m_[j];
		}
	}
	return out;
}

template<typename Type>
std::ostream& operator<<(std::ostream& flux, Vector<Type> const& vec){
	for(unsigned int i(0);i<vec.size();i++){
		flux<<vec(i)<<" ";
	}
	return flux;
}

template<typename Type>
std::istream& operator>>(std::istream& flux, Vector<Type>& vec){
	for(unsigned int i(0);i<vec.size();i++){
		flux>>vec(i);
	}
	return flux;
}

template<typename Type>
Vector<Type>& Vector<Type>::operator/=(Type const& d){
	for(unsigned int i(0);i<N_;i++){
		m_[i] /= d; 
	}
	return (*this);
}

template<typename Type>
Vector<Type> Vector<Type>::operator/(Type const& d) const{
	Vector<Type> tmp(*this);
	tmp /= d;
	return tmp;
}
/*}*/

/*methods that modify the class*/
/*{*/
template<>
inline Vector<double> Vector<double>::chop(double precision) const {
	Vector<double > tmp(*this);
	for(unsigned int i(0);i<N_;i++){
		if(std::abs(tmp.m_[i]) < precision ){tmp.m_[i]=0.0;}
	}
	return tmp;
}

template<>
inline Vector<std::complex<double> > Vector<std::complex<double> >::chop(double precision) const{
	Vector<std::complex<double> > tmp(*this);
	for(unsigned int i(0);i<N_;i++){
		if(std::abs(tmp.m_[i].imag()) < precision ){tmp.m_[i].imag(0.0);}
		if(std::abs(tmp.m_[i].real()) < precision ){tmp.m_[i].real(0.0);}
	}
	return tmp;
}

template<typename Type>
void Vector<Type>::set(){
	if(m_){ delete[] m_; }
	m_ = NULL;
	N_ = 0;
}

template<typename Type>
void Vector<Type>::set(unsigned int N){
	if(N_ != N){ 
		if(m_){ delete[] m_; }
		m_ = new Type[N];
		N_ = N;
	}
}

template<typename Type>
void Vector<Type>::set(unsigned int N, Type const& val){
	set(N);
	for(unsigned int i(0); i<N_; i++){
		m_[i] = val;
	}
}
/*}*/

template<typename Type>
Vector<Type> Vector<Type>::range(unsigned int min, unsigned int max) const {
	Vector<Type> out(max-min);
	for(unsigned int i(0);i<max-min;i++){
		out(i) = m_[min+i];
	}
	return out;
}

template<typename Type>
Type Vector<Type>::min() const {
	Type m(m_[0]);
	for(unsigned int i(1);i<N_;i++){
		if(m>m_[i]){ m=m_[i]; }
	}
	return m;
}

template<typename Type>
Type Vector<Type>::max() const {
	Type m(m_[0]);
	for(unsigned int i(1);i<N_;i++){
		if(m<m_[i]){ m=m_[i]; }
	}
	return m;
}

template<typename Type>
Type Vector<Type>::mean() const {
	double m(0);
	for(unsigned int i(0);i<N_;i++){
		m+=m_[i];
	}
	return m/N_;
}

/*Sort*/
/*{*/
template<typename Type>
bool Vector<Type>::is_sorted() const {
	for(unsigned int i(0);i<N_-1;i++) {
		if(m_[i]>m_[i+1]) {
			return true;
		}
	}
	return false;
}

template<typename Type>
void swap(Type& a, Type& b){
	Type tmp(a);
	a = b;
	b = tmp;
}

template<typename Type>
Vector<unsigned int> Vector<Type>::sort(){
	Vector<unsigned int> index(N_);
	for(unsigned int i(0);i<N_;i++) {
		index(i) = i;
	}
	while(is_sorted()) {
		for(unsigned int i(0);i<N_-1;i++) {
			if(m_[i]>m_[i+1]) {
				swap(m_[i],m_[i+1]);
				swap(index(i),index(i+1));
			}
		}
	}
	return index;
}

template<typename Type>
Vector<Type> Vector<Type>::sort(Vector<unsigned int> const& index) const{
	Vector<Type> out(index.size());
	for(unsigned int i(0);i<index.size();i++) {
		out(i) = m_[index(i)];
	}
	return out;
}
/*}*/
#endif
