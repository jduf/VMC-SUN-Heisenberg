#ifndef DEF_VECTOR_NEW
#define DEF_VECTOR_NEW

#include <iostream>
#include <cassert>
#include <complex>

/*!Class that implement a static array as a Vector
 *
 * - can be saved with Write.hpp 
 * - can be loaded with Read.hpp 
*/
template<typename Type>
class Vector : Matrix<Type>{
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

		/*!Accesses the ith entry of the Vector*/
		Type const& operator()(unsigned int const& i) const { assert(i<N_total); return m[i]; };
		/*!Sets the ith entry of the Vector*/
		Type& operator()(unsigned int const& i) { assert(i<N_total); return m[i]; };
};

/*constructors and destructor*/
/*{*/
template<typename Type>
Vector<Type>::Vector(): { }

template<typename Type>
Vector<Type>::Vector(unsigned int N): Matrix<Type>(N,1) { } 

template<typename Type>
Vector<Type>::Vector(unsigned int N, Type val): Matrix<Type>(N,1,val) { }

template<typename Type>
Vector<Type>::Vector(Vector<Type> const& vec): Matrix<Type>(vec) { }

template<typename Type>
Vector<Type>::~Vector(){ }
/*}*/
#endif
