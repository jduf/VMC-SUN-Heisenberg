#ifndef DEF_VECTOR
#define DEF_VECTOR

#include "Array.hpp"

#include <iostream>
#include <cassert>
#include <cmath> //allow abs(double) and abs(complex) 
#include <complex>

/*!Class that implement a static array as a vector
 *
 * - behaves as a std::vector<Type>
 * - can be saved with Write.hpp 
 * - can be loaded with Read.hpp 
 * */
template <typename Type>
class Vector : public Array<Type>{
	public:
		/*!Default constructor that initializes *v to NULL and N to 0*/
		Vector();
		/*!Initializes a static array of T of size N*/
		Vector(unsigned int N);
		/*!Initializes a static array of T of size N to a value val*/
		Vector(unsigned int N, Type val);
		/*!Deep copy*/
		Vector(Vector<Type> const& v);
		/*!Delete the static array*/
		~Vector();

		/*!Accesses the ith entry of the vector*/
		inline Type const& operator()(unsigned int const& i) const {assert(i<size()); return this->array[i];};
		/*!Set the ith entry of the vector*/
		inline Type& operator()(unsigned int const& i){assert(i<size()); return this->array[i];};
		/*!Multiplies this vector by a Type*/
		Vector<Type>& operator*=(Type const& d);
		Vector<Type> operator*(Type const& d) const;
		/*!Divide this vector by a Type*/
		Vector<Type>& operator/=(Type const& d);
		Vector<Type> operator/(Type const& d) const;
		/*!Scalar product*/
		Type operator*(Vector<Type> const& vec) const;
		/*!Exterior product*/
		//Matrice<Type> operator^(Vector<Type> const& vec) const;

		/*!Returns the norm of the vector*/
		double norm() const;

		/*!Returns the pointer on the array*/
		inline Type* ptr() const { return this->v; };
		/*!Returns the size of the vector*/
		inline unsigned int size() const { return this->total_size; };

		/*!Sets the entries to zero if they are close to 0*/
		void chop(double precision = 1e-10);

		/*!Set the all the values to val*/
		void set(Type val);

	private:
};

template <typename Type>
Vector<Type> operator*(Type const& d, Vector<Type> const& vec);
template <typename Type>
std::ostream& operator<<(std::ostream &flux, Vector<Type> const& v);

/*constructors and destructor*/
/*{*/
template <typename Type>
Vector<Type>::Vector(){ }

template <typename Type>
Vector<Type>::Vector(unsigned int N):
	Array<Type>(N)
{ }

template <typename Type>
Vector<Type>::Vector(unsigned int N, Type val):
	Array<Type>(N,val)
{ }

template <typename Type>
Vector<Type>::Vector(Vector<Type> const& vec):
	Array<Type>(vec)
{
}

template <typename Type>
Vector<Type>::~Vector(){ }
/*}*/

/*operators*/
/*{*/

template <typename Type>
Vector<Type>& Vector<Type>::operator*=(Type const& d){
	for(unsigned int i(0);i<this->total_size;i++){
		this->array[i] *= d;
	}
	return (*this);
}

template <typename Type>
Vector<Type> Vector<Type>::operator*(Type const& d) const{
	Vector<Type> vecout((*this));
	vecout *= d;
	return vecout;
}

template <typename Type>
Vector<Type> operator*(Type const& d, Vector<Type> const& vec){
	Vector<Type> vecout(vec);
	vecout *= d;	
	return vecout;
}

template <typename Type>
Type Vector<Type>::operator*(Vector<Type> const& vec) const{
	assert(this->total_size==vec.total_size);
	Type s(0.0);
	for(unsigned int i(0);i<this->total_size;i++){
		s += this->array[i]*vec.array[i];
	}
	return s;
}

template <typename Type>
Vector<Type>& Vector<Type>::operator/=(Type const& d){
	for(unsigned int i(0);i<this->total_size;i++){
		this->array[i] /= d;
	}
	return (*this);
}

template <typename Type>
Vector<Type> Vector<Type>::operator/(Type const& d) const{
	Vector<Type> vecout((*this));
	vecout /= d;	
	return vecout;
}

//template<typename Type>
//Matrice<Type> Vector<Type>::operator^(Vector<Type> const& vec) const{
	//assert(N == vec.N);
	//Matrice<Type> mat(N);
	//for(unsigned int i(0);i<N;i++){
		//for(unsigned int j(0);j<N;j++){
			//mat(i,j) = v[i]*vec.v[j];
		//}
	//}
	//return mat;
//}

template <typename Type>
std::ostream& operator<<(std::ostream &flux, Vector<Type> const& v){
	for(unsigned int i(0); i<v.size(); i++){
		flux << v(i) <<std::endl;
	}
	return flux;
}
/*}*/

/*methods that modify the class*/
/*{*/
template<>
inline void Vector<double>::chop(double precision){
	for(unsigned int i(0);i<this->total_size;i++){
		if(std::abs(this->array[i]) < precision ){this->array[i]=0.0;}
	}
}

template<>
inline void Vector<std::complex<double> >::chop(double precision){
	for(unsigned int i(0);i<this->total_size;i++){
		if(std::abs(this->array[i].real()) < precision ){this->array[i].real()=0.0;}
		if(std::abs(this->array[i].imag()) < precision ){this->array[i].imag()=0.0;}
	}
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template <typename Type>
double Vector<Type>::norm() const{
	Type d(0.0);
	for (unsigned int i(0);i<this->total_size; i++){
		d+= this->array[i]*this->array[i];
	}
	return d;
}
/*}*/
#endif
