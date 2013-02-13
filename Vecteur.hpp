#ifndef DEF_VECTEUR
#define DEF_VECTEUR

#include <iostream>
#include <cmath>
#include <cassert>

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

		/*!Sets the entries to zero if they are close to 0*/
		void chop();

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
		void fill_vecteur(T val);
};

template <typename T>
Vecteur<T> operator*(T const& d, Vecteur<T> const& vec);
template <typename T>
Vecteur<T> operator*(Vecteur<T> const& vec, T const& d);
template <typename T>
T operator*(Vecteur<T> const& vec1, Vecteur<T> const& vec2);
template <typename T>
std::ostream& operator<<(std::ostream &flux, Vecteur<T> const& v);

/*Constructors and destructor*/
/*{*/
template <typename T>
Vecteur<T>::Vecteur():v(NULL),N(0){}

template <typename T>
Vecteur<T>::Vecteur(unsigned int N):
	v(new T[N]),
	N(N)
{
	std::cout<<"taille : vecteur"<<std::endl;
}

template <typename T>
Vecteur<T>::Vecteur(unsigned int N, T val):
	v(new T[N]),
	N(N)
{
	std::cout<<"taille+const : vecteur"<<std::endl;
	fill_vecteur(val);
}

template <typename T>
Vecteur<T>::Vecteur(Vecteur<T> const& vec):
	v(new T[vec.size()]),
	N(vec.size())
{
	std::cout<<"copie : vecteur"<<std::endl;
	for(unsigned int i(0);i<N;i++){
			v[i] = vec(i);
	}
}

template <typename T>
Vecteur<T>::~Vecteur(){
	delete[] v;
	std::cout<<"destructeur : vecteur"<<std::endl;
}
/*}*/

/*operators*/
/*{*/
template <typename T>
Vecteur<T>& Vecteur<T>::operator=(Vecteur<T> const& vec){
	std::cout<<"affectation : vecteur"<<std::endl;
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
Vecteur<T>& Vecteur<T>::operator*=(T const& d){
	for(unsigned int i(0);i<N;i++){
		v[i] *= d;
	}
	return (*this);
}

template <typename T>
Vecteur<T> operator*(T const& d, Vecteur<T> const& vec){
	Vecteur<T> v(vec);
	v *= d;	
	return v;
}

template <typename T>
Vecteur<T> operator*(Vecteur<T> const& vec, T const& d){
	Vecteur<T> v(vec);
	v *= d;	
	return v;
}

template <typename T>
T operator*(Vecteur<T> const& vec1, Vecteur<T> const& vec2){
	assert(vec1.N==vec2.N);
	T s(0.0);
	for(unsigned int i(0);i<vec1.size();i++){
		s += vec1(i)*vec2(i);
	}
	return s;
}

template <typename T>
std::ostream& operator<<(std::ostream &flux, Vecteur<T> const& v){
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
void Vecteur<T>::fill_vecteur(T val){
	for(unsigned int i(0);i<N;i++){
		v[i] = val;
	}
}

template <typename T>
void Vecteur<T>::chop(){
	std::cerr<<"chop : attention utilise fabs... ne marche peut-être pas avec tous les types"<<std::endl;
	for(unsigned int i(0);i<N;i++){
		if(std::fabs(v[i]) < 1e-10 ){v[i]=0;}
	}
}
/*}*/

/*other methods*/
/*{*/
template <typename T>
void Vecteur<T>::print() const{
	for(unsigned int i(0);i<N;i++){
		std::cout<<v[i]<<std::endl;;
	}
}
/*}*/
#endif
