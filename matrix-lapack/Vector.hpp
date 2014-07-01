#ifndef DEF_VECTOR
#define DEF_VECTOR

#include "Matrix.hpp"
#include<functional>
#include<algorithm>

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
		Type const& operator()(unsigned int const& i) const 
		{assert(i<size_); return m_[i];}
		/*!Sets the (i,j)th entry of the Vector*/
		Type& operator()(unsigned int const& i) 
		{assert(i<size_); return m_[i];}

		/*!Deep copy assignment*/
		Vector<Type>& operator=(Vector<Type> vec); 
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
		/*!Multiplies a vectors by a scalar*/
		Vector<Type>& operator*=(Type const& d);
		Vector<Type> operator*(Type const& d) const;

		/*!Set the vector to 0*/
		void set();
		/*!Set the whole Vector to val*/
		void set(unsigned int N);
		/*!Set the vector to val*/
		void set(unsigned int N, Type const& val);

		/*!Sets the entries to zero if they are close to 0*/
		Vector<Type> chop(double precision = 1e-10) const;

		/*!Print the vector for vechevecica*/
		void print_mathematica() const;

		template<typename Function>
			void sort(Function cmp);
		template<typename Function>
			void sort(Function cmp, Vector<unsigned int>& index);
		template<typename Function>
			bool is_sorted(Function cmp) const;
		Vector<Type> order(Vector<unsigned int> const& index) const;
		Vector<Type> range(unsigned int min, unsigned int max) const;

		Type sum() const;
		Type mean() const;
		Type variance() const;
		Type max() const;
		Type min() const;

		/*!Returns the size of the Vector*/
		unsigned int size() const {return size_;}
		/*!Returns the pointer to the Vector*/
		Type* ptr() const {return m_;}

#ifdef DEF_IOFILES
		void header_rst(std::string const& s, RST& rst) const;
#endif

	protected:
		unsigned int size_; //!< number of rows
		Type *m_; //!< pointer to a static array

		//void set_null_pointer(){m_=NULL;}
		void swap_to_assign(Vector<Type>& v1,Vector<Type>& v2);
};

/*constructors and destructor*/
/*{*/
template<typename Type>
Vector<Type>::Vector():
	size_(0),
	m_(NULL)
{}

template<typename Type>
Vector<Type>::Vector(unsigned int N):
	size_(N),
	m_(size_?new Type[size_]:NULL)
{} 

template<typename Type>
Vector<Type>::Vector(unsigned int N, Type val):
	size_(N),
	m_(size_?new Type[size_]:NULL)
{
	for(unsigned int i(0); i<size_; i++){m_[i] = val;}
}

template<typename Type>
Vector<Type>::Vector(Vector<Type> const& vec):
	size_(vec.size_),
	m_(size_?new Type[size_]:NULL)
{
	for(unsigned int i(0);i<size_;i++){ m_[i] = vec.m_[i]; }
}

template<typename Type>
Vector<Type>::~Vector(){
	if(m_){ delete[]  m_; }
}

template<typename Type>
void Vector<Type>::swap_to_assign(Vector<Type>& v1,Vector<Type>& v2){
	std::swap(v1.m_,v2.m_);
	std::swap(v1.size_,v2.size_);
}
/*}*/

/*i/o methods*/
/*{*/
template<typename Type>
std::ostream& operator<<(std::ostream& flux, Vector<Type> const& v){
	for(unsigned int i(0);i<v.size();i++){
		flux<<v(i)<<" "; 
	}
	return flux;
}

template<typename Type>
std::istream& operator>>(std::istream& flux, Vector<Type> const& v){
	for(unsigned int i(0);i<v.size();i++){
		flux>>v.ptr()[i]; 
	}
	return flux;
}

#ifdef DEF_IOFILES
template<typename Type>
void Vector<Type>::header_rst(std::string const& s, RST& rst) const {
	rst.def(s,"Vector("+tostring(size_)+")"); 
}

template<typename Type>
IOFiles& operator<<(IOFiles& w, Vector<Type> const& v){
	if(w.is_binary()){
		w<<v.size();
		w.write(v.ptr(),v.size(),sizeof(Type));
	} else {
		w.stream()<<v;
	}
	return w;
}

template<typename Type>
IOFiles& operator>>(IOFiles& r, Vector<Type>& v){
	if(r.is_binary()){
		unsigned int size(0);
		r>>size;
		if(size != v.size()) {v.set(size);} 
		r.read(v.ptr(),v.size(),sizeof(Type));
	} else {
		r.stream()>>v;
	}
	return r;
}
#endif 
/*}*/

/*operators*/
/*{*/
template<typename Type>
Vector<Type>& Vector<Type>::operator=(Vector<Type> vec){
	swap_to_assign(*this,vec);
	return (*this);
}

template<typename Type>
Vector<Type>& Vector<Type>::operator+=(Vector<Type> const& vec){
	assert(size_ == vec.size_);
	for(unsigned int i(0);i<size_;i++){ m_[i] += vec.m_[i]; }
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
	assert(size_ == vec.size_);
	for(unsigned int i(0);i<size_;i++){ m_[i] -= vec.m_[i]; }
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
	for(unsigned int i(0);i<size_;i++){ out += m_[i] * vec.m_[i]; }
	return out;
}

template<typename Type>
Matrix<Type> Vector<Type>::operator^(Vector<Type> const& vec) const{
	Matrix<Type> out(vec.size(),vec.size());
	for(unsigned int i(0);i<size_;i++){
		for(unsigned int j(0);j<size_;j++){
			out(i,j) = m_[i] * vec.m_[j];
		}
	}
	return out;
}

template<typename Type>
Vector<Type>& Vector<Type>::operator/=(Type const& d){
	for(unsigned int i(0);i<size_;i++){ m_[i] /= d; }
	return (*this);
}

template<typename Type>
Vector<Type> Vector<Type>::operator/(Type const& d) const{
	Vector<Type> tmp(*this);
	tmp /= d;
	return tmp;
}

template<typename Type>
Vector<Type>& Vector<Type>::operator*=(Type const& d){
	for(unsigned int i(0);i<size_;i++){ m_[i] *= d; }
	return (*this);
}

template<typename Type>
Vector<Type> Vector<Type>::operator*(Type const& d) const{
	Vector<Type> tmp(*this);
	tmp *= d;
	return tmp;
}
/*}*/

/*methods that modify the class*/
/*{*/
template<>
inline Vector<double> Vector<double>::chop(double precision) const {
	Vector<double > tmp(*this);
	for(unsigned int i(0);i<size_;i++){
		if(std::abs(tmp.m_[i]) < precision ){tmp.m_[i]=0.0;}
	}
	return tmp;
}

template<>
inline Vector<std::complex<double> > Vector<std::complex<double> >::chop(double precision) const{
	Vector<std::complex<double> > tmp(*this);
	for(unsigned int i(0);i<size_;i++){
		if(std::abs(tmp.m_[i].imag()) < precision ){tmp.m_[i].imag(0.0);}
		if(std::abs(tmp.m_[i].real()) < precision ){tmp.m_[i].real(0.0);}
	}
	return tmp;
}

template<typename Type>
void Vector<Type>::set(){
	if(m_){ delete[] m_; }
	m_ = NULL;
	size_ = 0;
}

template<typename Type>
void Vector<Type>::set(unsigned int N){
	if(size_ != N){ 
		if(m_){ delete[] m_; }
		m_ = new Type[N];
		size_ = N;
	}
}

template<typename Type>
void Vector<Type>::set(unsigned int N, Type const& val){
	set(N);
	for(unsigned int i(0); i<size_; i++){
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

/*statistique*/
/*{*/
template<typename Type>
Type Vector<Type>::min() const {
	Type m(m_[0]);
	for(unsigned int i(1);i<size_;i++){
		if(m>m_[i]){ m=m_[i]; }
	}
	return m;
}

template<typename Type>
Type Vector<Type>::max() const {
	Type m(m_[0]);
	for(unsigned int i(1);i<size_;i++){
		if(m<m_[i]){ m=m_[i]; }
	}
	return m;
}

template<typename Type>
Type Vector<Type>::mean() const {
	double m(0.0);
	for(unsigned int i(0);i<size_;i++){ m+=m_[i]; }
	return m/size_;
}

template<typename Type>
Type Vector<Type>::variance() const {
	double m(mean());
	double v(0.0);
	for(unsigned int i(0);i<size_;i++){ v+=(m_[i]-m)*(m_[i]-m); }
	return v/size_;
}

template<typename Type>
Type Vector<Type>::sum() const {
	Type s(0.0);
	for(unsigned int i(0);i<size_;i++){ s += m_[i];}
	return s;
}
/*}*/

/*sort*/
/*{*/
template<typename Type>
template<typename Function>
void Vector<Type>::sort(Function cmp){
	std::sort(m_,m_+size_,cmp);
}

template<typename Type>
template<typename Function>
void Vector<Type>::sort(Function cmp, Vector<unsigned int>& index){
	index.set(size_);
	for(unsigned int i(0);i<size_;i++) { index(i) = i; }
	while(!is_sorted(cmp)) {
		for(unsigned int i(0);i<size_-1;i++) {
			if(!cmp(m_[i],m_[i+1])){
				std::swap(m_[i],m_[i+1]);
				std::swap(index(i),index(i+1));
			}
		}
	}
}

template<typename Type>
template<typename Function>
bool Vector<Type>::is_sorted(Function cmp) const {
	for(unsigned int i(0);i<size_-1;i++) {
		if(!cmp(m_[i],m_[i+1])){ return false; }
	}
	return true;
}

template<typename Type>
Vector<Type> Vector<Type>::order(Vector<unsigned int> const& index) const{
	Vector<Type> out(index.size());
	for(unsigned int i(0);i<index.size();i++) {
		out(i) = m_[index(i)];
	}
	return out;
}
/*}*/

/*double real(T)*/
/*{*/
inline double real(double const& x){ return x; }

inline double real(std::complex<double> const& x){ return std::real(x); }
/*}*/

/*double imag(T)*/
/*{*/
inline double imag(double const& x){ return x; }

inline double imag(std::complex<double> const& x){ return std::imag(x); }
/*}*/

/*double norm_squared(T)*/
/*{*/
inline double norm_squared(double x){ return x*x; }

inline double norm_squared(std::complex<double> x){ return std::norm(x); }
/*}*/
#endif
