#ifndef DEF_ARRAY
#define DEF_ARRAY

#include <iostream>
#include <cassert>
#include <complex>

template<typename Type>
class Array{
	public:
		/*!Default constructor that initializes *a to NULL and total_size to 0*/
		Array();
		/*!Initializes a static array of Type of size total_size*/
		Array(unsigned int total_size);
		/*!Initializes a static array of Type of size N*N*/
		Array(unsigned int total_size, Type val);
		/*!Deep copy*/
		Array(Array const& a);
		/*!Delete the static array*/
		~Array();

		/*!Does a deep copie*/
		Array& operator=(Array const& a);

		/*!Set the whole matrix to val*/
		void set(Type const& val);

	protected:
		Type *array; //!< pointer to a static array
		unsigned int total_size; //!< size of the array
};

/*constructors and destructor*/
/*{*/
template<typename Type>
Array<Type>::Array():
	array(NULL),
	total_size(0)
{ }

template<typename Type>
Array<Type>::Array(unsigned int total_size):
	array(new Type[total_size]),
	total_size(total_size)
{ } 

template<typename Type>
Array<Type>::Array(unsigned int total_size, Type val):
	array(new Type[total_size]),
	total_size(total_size)
{
	set(val);
}

template<typename Type>
Array<Type>::Array(Array<Type> const& a):
	array(new Type[a.total_size]),
	total_size(a.total_size)
{
	for(unsigned int i(0);i<total_size;i++){
		array[i] = a.array[i];
	}
}

template<typename Type>
Array<Type>::~Array(){
	if(array){
		delete[]  array;
		array = NULL;
	}
}
/*}*/

template <typename Type>
Array<Type>& Array<Type>::operator=(Array<Type> const& vec){
	if(this->total_size != vec.total_size){
		if(this->array){ delete[] this->array; }
		this->array = new Type[total_size];
		this->total_size = vec.total_size;
	}
	for(unsigned int i(0); i<total_size; i++){
		this->array[i] = vec.array[i];
	}
	return (*this);
}

template<typename Type>
void Array<Type>::set(Type const& val){
	for(unsigned int i(0); i<total_size; i++){
		array[i] = val;
	}
}
#endif
