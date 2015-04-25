#ifndef DEF_VECTOR
#define DEF_VECTOR

#include "Matrix.hpp"

/*{Description*/
/*!Class that implement a static array as a Vector
 *
 * - can be saved and loaded with IOFiles.hpp */
/*}*/
template<typename Type>
class Vector{
	public:
		/*!Default constructor that initializes *m to NULL and N to 0*/
		Vector();
		/*!Initializes a static array of Type of size N*/
		Vector(unsigned int N);
		/*!Initializes with {min,min+dx,...,max-dx,max}*/
		Vector(Type const& min, Type const& max, Type const& dx);
		/*!Initializes a static array of Type of size N to a value val*/
		Vector(unsigned int N, Type val);
		/*!Initialize with a std::vector*/
		Vector(std::vector<Type> const& vec);
		/*!Deep copy constructor*/
		Vector(Vector<Type> const& vec);
		/*!Move constructor (can't be default because of vec_)*/
		Vector(Vector<Type>&& vec);
		/*!Constructor that reads from file*/
		Vector(IOFiles& r);
		/*!Delete the static array*/
		~Vector();

		/*!Accesses the (i,j)th entry of the Vector*/
		Type const& operator()(unsigned int const& i) const 
		{ assert(i<size_); return vec_[i]; }
		/*!Sets the (i,j)th entry of the Vector*/
		Type& operator()(unsigned int const& i) 
		{ assert(i<size_); return vec_[i]; }

		/*!Assignment (using Copy-And-Swap Idiom)*/
		Vector<Type>& operator=(Vector<Type> vec); 
		/*!Additions this vector with another (m1 += m2 : m1 = m1+m2)*/
		Vector<Type>& operator+=(Vector<Type> const& vec);
		/*!Calls operator+=(Vector<Type> const& vec)*/
		Vector<Type> operator+(Vector<Type> const& vec) const;
		/*!Substracts this vector from another (m1 -= m2 : m1 = m1-m2)*/
		Vector<Type>& operator-=(Vector<Type> const& vec);
		/*!Calls operator-=(Vector<Type> const& vec)*/
		Vector<Type> operator-(Vector<Type> const& vec) const;
		/*!Multiplies two vectors (v1 *= v2 : v1 = v1*v2)*/
		Vector<Type> operator*(Vector<Type> const& vec) const;
		/*!Multiplies two vectors and get a matrix*/
		Matrix<Type> operator^(Vector<Type> const& vec) const;

		/*!Devides a vectors by a scalar*/
		Vector<Type>& operator/=(Type const& d);
		/*!Calls operator/=(Type const& d)*/
		Vector<Type> operator/(Type const& d) const;
		/*!Multiplies a vectors by a scalar*/
		Vector<Type>& operator*=(Type const& d);
		/*!Calls operator*=(Type const& d)*/
		Vector<Type> operator*(Type const& d) const;
		/*!Calls operator-=(Type const& d)*/
		Vector<Type>& operator-=(Type const& d);

		/*!Set the vector to 0*/
		void set();
		/*!Set the whole Vector to val*/
		void set(unsigned int N);
		/*!Set the vector to val*/
		void set(unsigned int N, Type const& val);

		/*!Sets the entries to zero if they are close to 0*/
		Vector<Type> chop(double precision = 1e-10) const;

		/*!Sort the datas according to std::sort(vec_,vec_+size_,cmp)*/
		template<typename Function>
			void sort(Function cmp);
		/*{Description*/
		/*!Sort the datas according cmp and my own implementation of sort. 
		 * \param cmp comparison method (std::greater<Type>(), std::less<Type>(),...)
		 * \param index stores the reordering */
		/*}*/
		template<typename Function>
			void sort(Function cmp, Vector<unsigned int>& index);
		/*!Returns true if the datas are sorted according cmp*/
		template<typename Function>
			bool is_sorted(Function cmp) const;
		/*!Returns a Vector composed by vec_ ordered by index*/
		Vector<Type> order(Vector<unsigned int> const& index) const;
		/*!Returns a vec_[min:max]*/
		Vector<Type> range(unsigned int min, unsigned int max) const;

		/*!Returns the sum over all the elements in vec_*/
		Type sum() const;
		/*!Returns the mean of the elements in vec_*/
		Type mean() const;
		/*!Returns the variance of the elements in vec_*/
		Type variance() const;
		/*!Returns the maximal value of vec_*/
		Type max() const;
		/*!Returns the minimal value of vec_*/
		Type min() const;
		/*!Returns the norm of the vector (sqrt(xx+yy+zz+...))*/
		double norm() const;

		/*!Returns the size of the Vector*/
		unsigned int size() const {return size_;}
		/*!Returns the pointer to the Vector*/
		Type* ptr() const {return vec_;}

		void header_rst(std::string const& s, RST& rst) const;

	private:
		unsigned int size_; //!< number of rows
		Type* vec_; //!< pointer to a static array

		/*!Copy-And-Swap Idiom*/
		void swap_to_assign(Vector<Type>& v1,Vector<Type>& v2);
};

/*constructors and destructor*/
/*{*/
template<typename Type>
Vector<Type>::Vector():
	size_(0),
	vec_(NULL)
{}

template<typename Type>
Vector<Type>::Vector(unsigned int N):
	size_(N),
	vec_(size_?new Type[size_]:NULL)
{} 

/*!Initializes with {min,min+dx,...,max-dx,max}*/
template<typename Type>
Vector<Type>::Vector(Type const& min, Type const& max, Type const& dx):
	size_(0),
	vec_(NULL)
{
	set((max-min)/dx+1);
	for(unsigned int i(0);i<size_;i++){ vec_[i] = min+i*dx; }
}

template<typename Type>
Vector<Type>::Vector(unsigned int N, Type val):
	size_(N),
	vec_(size_?new Type[size_]:NULL)
{
	for(unsigned int i(0);i<size_;i++){ vec_[i] = val; }
}

template<typename Type>
Vector<Type>::Vector(std::vector<Type> const& vec):
	size_(vec.size()),
	vec_(size_?new Type[size_]:NULL)
{
	std::copy(vec.begin(),vec.end(),vec_);
}

template<typename Type>
Vector<Type>::Vector(Vector<Type> const& vec):
	size_(vec.size_),
	vec_(size_?new Type[size_]:NULL)
{
	for(unsigned int i(0);i<size_;i++){ vec_[i] = vec.vec_[i]; }
	std::cout<<"vector copy"<<std::endl;
}

template<typename Type>
Vector<Type>::Vector(Vector<Type>&& vec):
	size_(vec.size_),
	vec_(vec.vec_)
{
	vec.vec_ = NULL;
	vec.size_ = 0;
}

template<typename Type>
Vector<Type>::Vector(IOFiles& r):
	size_(r.read<unsigned int>()),
	vec_(size_?new Type[size_]:NULL)
{
	r.read(vec_,size_,sizeof(Type));
}

template<typename Type>
Vector<Type>::~Vector(){
	if(vec_){ delete[]  vec_; }
}

template<typename Type>
void Vector<Type>::swap_to_assign(Vector<Type>& v1,Vector<Type>& v2){
	std::swap(v1.vec_,v2.vec_);
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
std::istream& operator>>(std::istream& flux, Vector<Type>& v){
	unsigned int size(v.size());
	for(unsigned int i(0);i<size;i++){
		flux>>v.ptr()[i]; 
	}
	return flux;
}

template<typename Type>
void Vector<Type>::header_rst(std::string const& s, RST& rst) const {
	rst.def(s,"Vector("+my::tostring(size_)+")"); 
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
	if(r.is_binary()){ v = std::move(Vector<Type>(r)); } 
	else { r.stream()>>v; }
	return r;
}
/*}*/

/*arithmetic operators*/
/*{*/
template<typename Type>
Vector<Type>& Vector<Type>::operator=(Vector<Type> vec){
	swap_to_assign(*this,vec);
	return (*this);
}

template<typename Type>
Vector<Type>& Vector<Type>::operator+=(Vector<Type> const& vec){
	assert(size_ == vec.size_);
	for(unsigned int i(0);i<size_;i++){ vec_[i] += vec.vec_[i]; }
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
	for(unsigned int i(0);i<size_;i++){ vec_[i] -= vec.vec_[i]; }
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
	for(unsigned int i(0);i<size_;i++){ out += vec_[i] * vec.vec_[i]; }
	return out;
}

template<typename Type>
Matrix<Type> Vector<Type>::operator^(Vector<Type> const& vec) const{
	Matrix<Type> out(vec.size(),vec.size());
	for(unsigned int i(0);i<size_;i++){
		for(unsigned int j(0);j<size_;j++){
			out(i,j) = vec_[i] * vec.vec_[j];
		}
	}
	return out;
}

template<typename Type>
Vector<Type>& Vector<Type>::operator/=(Type const& d){
	for(unsigned int i(0);i<size_;i++){ vec_[i] /= d; }
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
	for(unsigned int i(0);i<size_;i++){ vec_[i] *= d; }
	return (*this);
}

template<typename Type>
Vector<Type> Vector<Type>::operator*(Type const& d) const{
	Vector<Type> tmp(*this);
	tmp *= d;
	return tmp;
}

template<typename Type>
Vector<Type>& Vector<Type>::operator-=(Type const& d){
	for(unsigned int i(0);i<size_;i++){ vec_[i] -= d; }
	return (*this);
}
/*}*/

/*methods that modify the class*/
/*{*/
template<typename Type>
void Vector<Type>::set(){
	if(vec_){ delete[] vec_; }
	vec_ = NULL;
	size_ = 0;
}

template<typename Type>
void Vector<Type>::set(unsigned int N){
	if(size_ != N){
		if(vec_){ delete[] vec_; }
		vec_ = new Type[N];
		size_ = N;
	}
}

template<typename Type>
void Vector<Type>::set(unsigned int N, Type const& val){
	set(N);
	for(unsigned int i(0); i<size_; i++){ vec_[i] = val; }
}
/*}*/

/*methods that returns something*/
/*{*/
template<>
inline Vector<double> Vector<double>::chop(double precision) const {
	Vector<double > tmp(*this);
	for(unsigned int i(0);i<size_;i++){
		if(std::abs(tmp.vec_[i]) < precision ){tmp.vec_[i]=0.0;}
	}
	return tmp;
}

template<>
inline Vector<std::complex<double> > Vector<std::complex<double> >::chop(double precision) const{
	Vector<std::complex<double> > tmp(*this);
	for(unsigned int i(0);i<size_;i++){
		if(std::abs(tmp.vec_[i].real()) < precision ){tmp.vec_[i].real(0.0);}
		if(std::abs(tmp.vec_[i].imag()) < precision ){tmp.vec_[i].imag(0.0);}
	}
	return tmp;
}

template<typename Type>
Vector<Type> Vector<Type>::range(unsigned int min, unsigned int max) const {
	Vector<Type> out(max-min);
	for(unsigned int i(0);i<max-min;i++){
		out.vec_[i] = vec_[min+i];
	}
	return out;
}
/*}*/

/*statistique*/
/*{*/
template<typename Type>
Type Vector<Type>::min() const {
	Type m(vec_[0]);
	for(unsigned int i(1);i<size_;i++){
		if(m>vec_[i]){ m=vec_[i]; }
	}
	return m;
}

template<typename Type>
Type Vector<Type>::max() const {
	Type m(vec_[0]);
	for(unsigned int i(1);i<size_;i++){
		if(m<vec_[i]){ m=vec_[i]; }
	}
	return m;
}

template<typename Type>
Type Vector<Type>::mean() const {
	double m(0.0);
	for(unsigned int i(0);i<size_;i++){ m+=vec_[i]; }
	return m/size_;
}

template<typename Type>
Type Vector<Type>::variance() const {
	double m(mean());
	double v(0.0);
	for(unsigned int i(0);i<size_;i++){ v+=(vec_[i]-m)*(vec_[i]-m); }
	return v/size_;
}

template<typename Type>
Type Vector<Type>::sum() const {
	Type s(0.0);
	for(unsigned int i(0);i<size_;i++){ s += vec_[i];}
	return s;
}

template<typename Type>
double Vector<Type>::norm() const {
	double norm(0);
	for(unsigned int i(0);i<size_;i++){ norm += norm_squared(vec_[i]);}
	return sqrt(norm);
}
/*}*/

/*sort*/
/*{*/
template<typename Type>
template<typename Function>
void Vector<Type>::sort(Function cmp){
	std::sort(vec_,vec_+size_,cmp);
}

template<typename Type>
template<typename Function>
void Vector<Type>::sort(Function cmp, Vector<unsigned int>& index){
	index.set(size_);
	for(unsigned int i(0);i<size_;i++){ index(i) = i; }
	while(!is_sorted(cmp)) {
		for(unsigned int i(0);i<size_-1;i++){
			if(!cmp(vec_[i],vec_[i+1])){
				std::swap(vec_[i],vec_[i+1]);
				std::swap(index(i),index(i+1));
			}
		}
	}
}

template<typename Type>
template<typename Function>
bool Vector<Type>::is_sorted(Function cmp) const {
	for(unsigned int i(0);i<size_-1;i++) {
		if(!cmp(vec_[i],vec_[i+1])){ return false; }
	}
	return true;
}

template<typename Type>
Vector<Type> Vector<Type>::order(Vector<unsigned int> const& index) const{
	assert(size_ == index.size());
	Vector<Type> out(size_);
	for(unsigned int i(0);i<size_;i++){ out(i) = vec_[index(i)]; }
	return out;
}
/*}*/

namespace my {
	template<typename Type>
		bool are_equal(Vector<Type> const& x, Vector<Type> const& y, double abs_tol=1e-14, double rel_tol=1e-14){
			if(x.size() != y.size()){ return false; }
			else {
				for(unsigned int i(0);i<x.size();i++){
					if(!are_equal(x(i),y(i),abs_tol,rel_tol)){ return false ; }
				}
				return true;
			}
		}
}
#endif
