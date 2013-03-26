#ifndef DEF_VECTEUR
#define DEF_VECTEUR

#include <iostream>
#include <cassert>

template <typename T>
class Matrice;

/*!Class that implement a static array as a vector
 *
 * - behaves as a std::vector<T>
 * - can be saved with Write.hpp 
 * - can be loaded with Read.hpp 
 * */
template <typename T>
class Vecteur{
	public:
		/*!Default constructor that initializes *v to NULL and N to 0*/
		Vecteur();
		/*!Initializes a static array of T of size N*/
		Vecteur(unsigned int N);
		/*!Initializes a static array of T of size N to a value val*/
		Vecteur(unsigned int N, T val);
		/*!Deep copy*/
		Vecteur(Vecteur<T> const& v);
		/*!Delete the static array*/
		~Vecteur();

		/*!Does a deep copie*/
		Vecteur& operator=(Vecteur const& mat);
		/*!Accesses the ith entry of the vector*/
		inline T const& operator()(unsigned int const& i) const {assert(i<N); return v[i];};
		/*!Set the ith entry of the vector*/
		inline T& operator()(unsigned int const& i){assert(i<N); return v[i];};
		/*!Multiplies this vector by a double*/
		Vecteur<T>& operator*=(T const& d);
		Vecteur<T> operator*(T const& d) const;
		/*!Scalar product*/
		T operator*(Vecteur<T> const& vec) const;
		/*!Exterior product*/
		Matrice<T> operator^(Vecteur<T> const& vec) const;

		/*!Returns the pointer on the array*/
		inline T* ptr() const { return v; };
		/*!Returns the size of the vector*/
		inline unsigned int size() const { return N; };

	private:
		T *v; //!< pointer to a static array
		unsigned int N; //!< size of the static array
		
		/*!Set the all the values to val*/
		void fill_vecteur(T val);
};

template <typename T>
Vecteur<T> operator*(T const& d, Vecteur<T> const& vec);
template <typename T>
std::ostream& operator<<(std::ostream &flux, Vecteur<T> const& v);

/*Constructors and destructor*/
/*{*/
template <typename T>
Vecteur<T>::Vecteur():v(NULL),N(0){ }

template <typename T>
Vecteur<T>::Vecteur(unsigned int N):
	v(new T[N]),
	N(N)
{ }

template <typename T>
Vecteur<T>::Vecteur(unsigned int N, T val):
	v(new T[N]),
	N(N)
{
	fill_vecteur(val);
}

template <typename T>
Vecteur<T>::Vecteur(Vecteur<T> const& vec):
	v(new T[vec.size()]),
	N(vec.size())
{
	for(unsigned int i(0);i<N;i++){
			v[i] = vec(i);
	}
}

template <typename T>
Vecteur<T>::~Vecteur(){
	if(v){
		delete[] v;
		v = NULL;
	}
}
/*}*/

/*operators*/
/*{*/
template <typename T>
Vecteur<T>& Vecteur<T>::operator=(Vecteur<T> const& vec){
	if(this->N != vec.N){
		if(this->v){ delete[] this->v; }
		this->v = new T[N];
		this->N = vec.N;
	}
	for(unsigned int i(0); i<N; i++){
		this->v[i] = vec.v[i];
	}
	return (*this);
}

template <typename T>
Vecteur<T>& Vecteur<T>::operator*=(T const& d){
	for(unsigned int i(0);i<N;i++){
		v[i] *= d;
	}
	return (*this);
}

template <typename T>
Vecteur<T> Vecteur<T>::operator*(T const& d) const{
	Vecteur<T> vecout((*this));
	vecout *= d;	
	return vecout;
}

template <typename T>
Vecteur<T> operator*(T const& d, Vecteur<T> const& vec){
	Vecteur<T> vecout(vec);
	vecout *= d;	
	return vecout;
}

template <typename T>
T Vecteur<T>::operator*(Vecteur<T> const& vec) const{
	assert(N==vec.N);
	T s(0.0);
	for(unsigned int i(0);i<N;i++){
		s += v[i]*vec.v[i];
	}
	return s;
}

template<typename T>
Matrice<T> Vecteur<T>::operator^(Vecteur<T> const& vec) const{
	assert(N == vec.N);
	Matrice<T> mat(N);
	for(unsigned int i(0);i<N;i++){
		for(unsigned int j(0);j<N;j++){
			mat(i,j) = v[i]*vec.v[j];
		}
	}
	return mat;
}

template <typename T>
std::ostream& operator<<(std::ostream &flux, Vecteur<T> const& v){
	for(unsigned int i(0); i<v.size(); i++){
		flux << v(i) <<std::endl;
	}
	return flux;
}
/*}*/

/*methods that modify the class*/
/*{*/
template <typename T>
void Vecteur<T>::fill_vecteur(T val){
	for(unsigned int i(0);i<N;i++){
		v[i] = val;
	}
}
/*}*/
#endif
