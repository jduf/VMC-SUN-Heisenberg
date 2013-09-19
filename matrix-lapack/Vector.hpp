#ifndef DEF_Vector
#define DEF_Vector

#include <iostream>
#include <cassert>
#include <cvech> //allow abs(double) and abs(complex) 
#include <complex>

/*!Class that implement a static array as a Vector
 *
 * - can be saved with Write.hpp 
 * - can be loaded with Read.hpp 
*/
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
		Type const& operator()(unsigned int const& i) const { assert(i<N); return m[i]; };
		/*!Sets the (i,j)th entry of the Vector*/
		Type& operator()(unsigned int const& i) { assert(i<N); return m[i]; };

		/*!Deep copy assignment*/
		Vector<Type>& operator=(Vector<Type> const& vec); 
		/*!Additions this vecrice with another*/
		Vector<Type>& operator+=(Vector<Type> const& vec);
		Vector<Type> operator+(Vector<Type> const& vec) const;
		/*!Substracts this vecrice from another (m1 -= m2 : m1 = m1-m2)*/
		Vector<Type>& operator-=(Vector<Type> const& vec);
		Vector<Type> operator-(Vector<Type> const& vec) const;
		/*!Multiplies two vecrices (m1 *= m2 : m1 = m1*m2)*/
		Vector<Type>& operator*=(Vector<Type> const& vec);
		Vector<Type> operator*(Vector<Type> const& vec) const;

		/*!Set the whole Vector to val*/
		void set(Type const& val);
		/*!Sets the entries to zero if they are close to 0*/
		void chop(double precision = 1e-10);
		/*!Print the vector for vechevecica*/
		void print_mathematica();

		/*!Returns the transpose of any Vector*/
		Vector<Type> transpose() const;
		/*!Returns the conjugate transpose of complex Vector (may give an error) */
		Vector<Type> trans_conj() const;
		/*!Returns the diagonal elements in an vector*/
		Vector<Type> diag() const;

		/*!Returns the pointer to the Vector*/
		Type* ptr() const { return m; }
		/*!Returns the size of the Vector*/
		unsigned int size() const { return N; }

	protected:
		Type *m; //!< pointer to a static array
		unsigned int N; //!< number of rows

		void set_null_pointer(){m=NULL;}
};

template<typename Type>
std::ostream& operator<<(std::ostream& flux, Vector<Type> const& vec);
template<typename Type>
std::istream& operator>>(std::istream& flux, Vector<Type>& vec);

/*constructors and destructor*/
/*{*/
template<typename Type>
Vector<Type>::Vector():
	m(NULL),
	N(0)
{
	//std::cout<<"called default"<<std::endl;
}

template<typename Type>
Vector<Type>::Vector(unsigned int N):
	m(new Type[N]),
	N(N)
{
	//std::cout<<"called N,N_col"<<std::endl;
} 

template<typename Type>
Vector<Type>::Vector(unsigned int N, Type val):
	m(new Type[N]),
	N(N)
{
	//std::cout<<"called N,N_col,val"<<std::endl;
	set(val);
}

template<typename Type>
Vector<Type>::Vector(Vector<Type> const& vec):
	m(new Type[vec.N]),
	N(vec.N)
{
	for(unsigned int i(0);i<N;i++){
		this->m[i] = vec.m[i];
	}
}

template<typename Type>
Vector<Type>::~Vector(){
	if(m){
		//std::cout<<" with delete";
		delete[]  m;
		m = NULL;
	}
}
/*}*/

/*operators*/
/*{*/
template<typename Type>
Vector<Type>& Vector<Type>::operator=(Vector<Type> const& vec){
	if(this->N_col != vec.N_col ||  this->N != vec.N){
		if(this->m){ delete[] this->m;}
		this->m = new Type[vec.N*vec.N_col];
		this->N_total = vec.N_total;
		this->N = vec.N;
		this->N_col = vec.N_col;
	}
	for(unsigned int i(0); i<this->N_total; i++){
		this->m[i] = vec.m[i];
	}
	return (*this);
}

template<typename Type>
Vector<Type>& Vector<Type>::operator+=(Vector<Type> const& vec){
	assert(N == vec.N && N_col == vec.N_col);
	for(unsigned int i(0);i<N_total;i++){
		this->m[i] += vec.m[i];
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
	assert(N == vec.N && N_col == vec.N_col);
	for(unsigned int i(0);i<N_total;i++){
		this->m[i] -= vec.m[i];
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
Vector<Type>& Vector<Type>::operator*=(Vector<Type> const& vec){
	assert(this->N_col==vec.N);
	Vector<Type> tmp(*this);
	if(this->N_col!=this->N || vec.N_col!=vec.N ){
		if(this->m){ delete[] this->m;}
		this->m = new Type[this->N*vec.N_col];
		this->N_total = this->N*vec.N_col;
		this->N_col = vec.N_col;
	}
	for(unsigned int i(0);i<this->N;i++){
		for(unsigned int j(0);j<this->N_col;j++){
			this->m[i+j*this->N] = 0.0;
			for(unsigned int k(0);k<tmp.N_col;k++){
				//this->m[i+j*N] += tmp(i,k) * vec(k,j);
				this->m[i+j*this->N] += tmp.m[i+k*tmp.N] * vec.m[k+j*vec.N];
			}
		}
	}
	return (*this);
}

template<typename Type>
Vector<Type> Vector<Type>::operator*(Vector<Type> const& vec) const{
	assert(this->N_col==vec.N);
	Vector<Type> vecout(this->N,vec.N_col);
	for(unsigned int i(0);i<vecout.N;i++){
		for(unsigned int j(0);j<vecout.N_col;j++){
			vecout.m[i+j*vecout.N] = 0.0;
			for(unsigned int k(0);k<vec.N;k++){
				//vecout.m[i+j*vecout.N] += (*this)(i,k) * vec(k,j);
				vecout.m[i+j*vecout.N] += this->m[i+k*this->N] * vec.m[k+j*vec.N];
			}
		}
	}
	return vecout;
}

template<typename Type>
std::ostream& operator<<(std::ostream& flux, Vector<Type> const& vec){
	for(unsigned int i(0);i<vec.row();i++){
		for(unsigned int j(0);j<vec.col();j++){
			flux<<vec(i,j)<<" ";
		}
		flux<<std::endl;
	}
	return flux;
}

template<typename Type>
std::istream& operator>>(std::istream& flux, Vector<Type>& vec){
	for(unsigned int i(0);i<vec.row();i++){
		for(unsigned int j(0);j<vec.col();j++){
			flux>>vec(i,j);
		}
	}
	return flux;
}
/*}*/

/*methods that modify the class*/
/*{*/
template<>
inline void Vector<double>::chop(double precision){
	for(unsigned int i(0);i<N_total;i++){
		if(std::abs(m[i]) < precision ){m[i]=0.0;}
	}
}

template<>
inline void Vector<std::complex<double> >::chop(double precision){
	for(unsigned int i(0);i<N_total;i++){
		if(std::abs(m[i].imag()) < precision ){m[i].imag()=0.0;}
		if(std::abs(m[i].real()) < precision ){m[i].real()=0.0;}
	}
}

template<typename Type>
void Vector<Type>::set(Type const& val){
	for(unsigned int i(0); i<N_total; i++){
		m[i] = val;
	}
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename Type>
Vector<Type> Vector<Type>::transpose() const{
	Vector<Type> tmp(this->N_col,this->N);
	for(unsigned int i(0);i<N_col;i++){
		for(unsigned int j(0);j<N;j++){
			tmp.m[i+j*tmp.N] = m[j+i*this->N];
		}
	}
	return tmp;
}

template<typename Type>
Vector<Type> Vector<Type>::trans_conj() const{
	Vector<Type> tmp(this->N_col,this->N);
	for(unsigned int i(0);i<N_col;i++){
		for(unsigned int j(0);j<N;j++){
			tmp.m[i+j*tmp.N] = conj(m[j+i*this->N]);
		}
	}
	return tmp;
}

template<typename Type>
Vector<Type> Vector<Type>::diag() const{
	unsigned int N(0);
	if(this->N < this->N_col){N=N_col;}
	else{N=N;}
	Vector<Type> v(N,1);
	for(unsigned int i(0);i<N;i++){
		v[i] = this->m[i*(this->N+1)];
	}
	return v;
}
/*}*/
#endif
