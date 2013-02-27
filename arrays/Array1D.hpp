#ifndef DEF_ARRAY1D
#define DEF_ARRAY1D

#include <iostream>
#include <cassert>

template <typename T>
class Array1D{
	public:
		/*!Default constructor that initializes *v to NULL and N to 0*/
		Array1D();
		/*!Initializes a static array of T of size N*/
		Array1D(unsigned int N);
		/*!Initializes a static array of T of size N to a value val*/
		Array1D(unsigned int N, T val);
		/*!Deep copy*/
		Array1D(Array1D<T> const& a);
		/*!Delete the static array*/
		~Array1D();

		/*!Does a deep copie*/
		Array1D& operator=(Array1D const& mat);
		/*!Accesses the ith entry of the vector*/
		inline T const& operator()(unsigned int const& i) const {assert(i<N); return v[i];};
		/*!Set the ith entry of the vector*/
		inline T& operator()(unsigned int const& i){assert(i<N); return v[i];};

		/*!Returns the pointer on the array*/
		inline T* ptr() const { return v; };
		/*!Returns the size of the vector*/
		inline unsigned int size() const { return N; };

		void test() const;
		void print() const;

	private:
		T *v; //!< pointer to a static array
		unsigned int N; //!< size of the static array
		
		/*!Set the all the values to val*/
		void fill_Array1D(T val);
};

template <typename T>
std::ostream& operator<<(std::ostream &flux, Array1D<T> const& v);

/*Constructors and destructor*/
/*{*/
template <typename T>
Array1D<T>::Array1D():v(NULL),N(0){}

template <typename T>
Array1D<T>::Array1D(unsigned int N):
	v(new T[N]),
	N(N)
{
	std::cout<<"taille : Array1D"<<std::endl;
}

template <typename T>
Array1D<T>::Array1D(unsigned int N, T val):
	v(new T[N]),
	N(N)
{
	std::cout<<"taille+const : Array1D"<<std::endl;
	fill_Array1D(val);
}

template <typename T>
Array1D<T>::Array1D(Array1D<T> const& vec):
	v(new T[vec.size()]),
	N(vec.size())
{
	std::cout<<"copie : Array1D"<<std::endl;
	for(unsigned int i(0);i<N;i++){
			v[i] = vec(i);
	}
}

template <typename T>
Array1D<T>::~Array1D(){
	delete[] v;
	std::cout<<"destructeur : Array1D"<<std::endl;
}
/*}*/

/*operators*/
/*{*/
template <typename T>
Array1D<T>& Array1D<T>::operator=(Array1D<T> const& vec){
	std::cout<<"affectation : Array1D"<<std::endl;
	if(this->N != vec.N){
		delete[] this->v;
		this->N = vec.N;
	}
	for(unsigned int i(0); i<N; i++){
		this->v[i] = vec.v[i];
	}
	return (*this);
}

template <typename T>
std::ostream& operator<<(std::ostream &flux, Array1D<T> const& v){
	for(unsigned int i(0); i<v.size()-1; i++){
		flux << v(i) <<std::endl;
	}
	flux << v(v.size()-1);
	return flux;
}
/*}*/

/*methods that modify the class*/
/*{*/
template <typename T>
void Array1D<T>::fill_Array1D(T val){
	for(unsigned int i(0);i<N;i++){
		v[i] = val;
	}
}
/*}*/

/*other methods*/
/*{*/
template <typename T>
void Array1D<T>::print() const{
	for(unsigned int i(0);i<N;i++){
		std::cout<<v[i]<<std::endl;;
	}
}
/*}*/
#endif
