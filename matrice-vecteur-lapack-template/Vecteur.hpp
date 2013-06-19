#ifndef DEF_VECTEUR
#define DEF_VECTEUR

#include <iostream>
#include <cassert>
#include <cmath> //allow abs(double) and abs(complex) 
#include <complex>

template <typename Type>
class Matrice;

/*!Class that implement a static array as a vector
 *
 * - behaves as a std::vector<Type>
 * - can be saved with Write.hpp 
 * - can be loaded with Read.hpp 
 * */
template <typename Type>
class Vecteur{
	public:
		/*!Default constructor that initializes *v to NULL and N to 0*/
		Vecteur();
		/*!Initializes a static array of T of size N*/
		Vecteur(unsigned int N);
		/*!Initializes a static array of T of size N to a value val*/
		Vecteur(unsigned int N, Type val);
		/*!Deep copy*/
		Vecteur(Vecteur<Type> const& v);
		/*!Delete the static array*/
		~Vecteur();

		/*!Does a deep copie*/
		Vecteur& operator=(Vecteur const& vec);
		/*!Accesses the ith entry of the vector*/
		inline Type const& operator()(unsigned int const& i) const {assert(i<N); return v[i];};
		/*!Set the ith entry of the vector*/
		inline Type& operator()(unsigned int const& i){assert(i<N); return v[i];};
		/*!Multiplies this vector by a Type*/
		Vecteur<Type>& operator*=(Type const& d);
		Vecteur<Type> operator*(Type const& d) const;
		/*!Divide this vector by a Type*/
		Vecteur<Type>& operator/=(Type const& d);
		Vecteur<Type> operator/(Type const& d) const;
		/*!Scalar product*/
		Type operator*(Vecteur<Type> const& vec) const;
		/*!Exterior product*/
		Matrice<Type> operator^(Vecteur<Type> const& vec) const;

		/*!Returns the norm of the vector*/
		double norm() const;

		/*!Returns the pointer on the array*/
		inline Type* ptr() const { return v; };
		/*!Returns the size of the vector*/
		inline unsigned int size() const { return N; };

		/*!Sets the entries to zero if they are close to 0*/
		void chop(double precision = 1e-10);

		/*!Set the all the values to val*/
		void set(Type val);

	private:
		Type *v; //!< pointer to a static array
		unsigned int N; //!< size of the static array
		
};

template <typename Type>
Vecteur<Type> operator*(Type const& d, Vecteur<Type> const& vec);
template <typename Type>
std::ostream& operator<<(std::ostream &flux, Vecteur<Type> const& v);

/*constructors and destructor*/
/*{*/
template <typename Type>
Vecteur<Type>::Vecteur():v(NULL),N(0){ }

template <typename Type>
Vecteur<Type>::Vecteur(unsigned int N):
	v(new Type[N]),
	N(N)
{ }

template <typename Type>
Vecteur<Type>::Vecteur(unsigned int N, Type val):
	v(new Type[N]),
	N(N)
{
	set(val);
}

template <typename Type>
Vecteur<Type>::Vecteur(Vecteur<Type> const& vec):
	v(new Type[vec.size()]),
	N(vec.size())
{
	for(unsigned int i(0);i<N;i++){
			v[i] = vec(i);
	}
}

template <typename Type>
Vecteur<Type>::~Vecteur(){
	if(v){
		delete[] v;
		v = NULL;
	}
}
/*}*/

/*operators*/
/*{*/
template <typename Type>
Vecteur<Type>& Vecteur<Type>::operator=(Vecteur<Type> const& vec){
	if(this->N != vec.N){
		if(this->v){ delete[] this->v; }
		this->v = new Type[N];
		this->N = vec.N;
	}
	for(unsigned int i(0); i<N; i++){
		this->v[i] = vec.v[i];
	}
	return (*this);
}

template <typename Type>
Vecteur<Type>& Vecteur<Type>::operator*=(Type const& d){
	for(unsigned int i(0);i<N;i++){
		v[i] *= d;
	}
	return (*this);
}

template <typename Type>
Vecteur<Type> Vecteur<Type>::operator*(Type const& d) const{
	Vecteur<Type> vecout((*this));
	vecout *= d;	
	return vecout;
}

template <typename Type>
Vecteur<Type> operator*(Type const& d, Vecteur<Type> const& vec){
	Vecteur<Type> vecout(vec);
	vecout *= d;	
	return vecout;
}

template <typename Type>
Type Vecteur<Type>::operator*(Vecteur<Type> const& vec) const{
	assert(N==vec.N);
	Type s(0.0);
	for(unsigned int i(0);i<N;i++){
		s += v[i]*vec.v[i];
	}
	return s;
}

template <typename Type>
Vecteur<Type>& Vecteur<Type>::operator/=(Type const& d){
	for(unsigned int i(0);i<N;i++){
		v[i] /= d;
	}
	return (*this);
}

template <typename Type>
Vecteur<Type> Vecteur<Type>::operator/(Type const& d) const{
	Vecteur<Type> vecout((*this));
	vecout /= d;	
	return vecout;
}

template<typename Type>
Matrice<Type> Vecteur<Type>::operator^(Vecteur<Type> const& vec) const{
	assert(N == vec.N);
	Matrice<Type> mat(N);
	for(unsigned int i(0);i<N;i++){
		for(unsigned int j(0);j<N;j++){
			mat(i,j) = v[i]*vec.v[j];
		}
	}
	return mat;
}

template <typename Type>
std::ostream& operator<<(std::ostream &flux, Vecteur<Type> const& v){
	for(unsigned int i(0); i<v.size(); i++){
		flux << v(i) <<std::endl;
	}
	return flux;
}
/*}*/

/*methods that modify the class*/
/*{*/
template <typename Type>
void Vecteur<Type>::set(Type val){
	for(unsigned int i(0);i<N;i++){
		v[i] = val;
	}
}

template<>
inline void Vecteur<double>::chop(double precision){
	for(unsigned int i(0);i<N;i++){
		if(std::abs(v[i]) < precision ){v[i]=0.0;}
	}
}

template<>
inline void Vecteur<std::complex<double> >::chop(double precision){
	for(unsigned int i(0);i<N*N;i++){
		if(std::abs(v[i].real()) < precision ){v[i].real()=0.0;}
		if(std::abs(v[i].imag()) < precision ){v[i].imag()=0.0;}
	}
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template <typename Type>
double Vecteur<Type>::norm() const{
	Type d(0.0);
	for (unsigned int i(0);i<N;i++){
		d+= v[i]*v[i];
	}
	return d;
}
/*}*/
#endif
