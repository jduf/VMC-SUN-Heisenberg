#ifndef DEF_VECTEUR
#define DEF_VECTEUR

#include <iostream>
#include <cmath>

template <typename T>
class Vecteur{
	public:
/*Constructors and destructor*/
		Vecteur();
		Vecteur(unsigned int N);
		Vecteur(unsigned int N, T val);
		//Vecteur(unsigned int N, V val);
		Vecteur(Vecteur<T> const& v);
		~Vecteur();

/*operators*/
		Vecteur& operator=(Vecteur const& mat);
		inline T const& operator()(unsigned int const& i) const {return v[i];};
		inline T& operator()(unsigned int const& i){return v[i];};
		Vecteur<T>& operator*=(T const& d);

/*methods that modify the class*/
		void chop();

/*other methods*/
		inline T* ptr() const { return v; };
		inline unsigned int size() const { return N; };

		void test() const;
		void print() const;

	private:
		T *v;
		unsigned int N;
		
/*methods that modify the class*/
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
{ }

template <typename T>
Vecteur<T>::Vecteur(unsigned int N, T val):
	v(new T[N]),
	N(N)
{ fill_vecteur(val); }

template <typename T>
Vecteur<T>::Vecteur(Vecteur const& vec):
	v(new T[vec.size()]),
	N(vec.size())
{
	for(unsigned int i(0);i<N;i++){
			v[i] = vec(i);
	}
}

template <typename T>
Vecteur<T>::~Vecteur(){
	delete[] v;
}
/*}*/

/*operators*/
/*{*/
template <typename T>
Vecteur<T>& Vecteur<T>::operator=(Vecteur<T> const& vec){
	if(this->N!=vec.N){
		std::cerr<<"impossible d'affecter deux Vecteurs si dim1 =! dim2"<<std::endl;
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
	if(vec1.size() != vec2.size()){
		std::cerr<<"scalar product between vector of different size"<<std::endl;
	}
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
