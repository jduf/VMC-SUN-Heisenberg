#ifndef DEF_MATRIX
#define DEF_MATRIX

template<typename Type>
class Vector;

#include <cassert>
#include "IOFiles.hpp"

/*{*/
/*!Class that implement a static array as a Matrix
 *
 * - can be saved and loaded with IOFiles.hpp */
/*}*/
template<typename Type>
class Matrix{
	public:
		/*!Default constructor that initializes *m to NULL and N to 0*/
		Matrix() = default;
		/*!Initializes a static array of Type of size N_row*N_col*/
		Matrix(unsigned int N_row, unsigned int N_col);
		/*!Initializes a static array of Type of size N_row*N_col to a value val*/
		Matrix(unsigned int N_row, unsigned int N_col, Type val);
		/*!Deep copy constructor*/
		Matrix(Matrix<Type> const& mat);
		/*!Move constructor (can't be default because of mat_)*/
		Matrix(Matrix<Type>&& mat);
		/*!Constructor that reads from file*/
		Matrix<Type>(IOFiles& r);
		/*!Delete the static array*/
		virtual ~Matrix();

		/*!Accesses the (i,j)th entry of the matrix*/
		Type const& operator()(unsigned int const& i, unsigned int const& j) const
		{ assert(i<row_ && j<col_); return mat_[i+j*row_]; }
		/*!Sets the (i,j)th entry of the matrix*/
		Type& operator()(unsigned int const& i, unsigned int const& j)
		{ assert(i<row_ && j<col_); return mat_[i+j*row_]; }

		/*!Assignment (using Copy-And-Swap Idiom)*/
		Matrix<Type>& operator=(Matrix<Type> mat);
		/*!Additions this matrice with another*/
		Matrix<Type>& operator+=(Matrix<Type> const& mat);
		/*!Calls operator+=(Matrix<Type> const& mat)*/
		Matrix<Type> operator+(Matrix<Type> const& mat) const;
		/*!Substracts this matrice from another (m1 -= m2 : m1 = m1-m2)*/
		Matrix<Type>& operator-=(Matrix<Type> const& mat);
		/*!Calls operator-=(Matrix<Type> const& mat)*/
		Matrix<Type> operator-(Matrix<Type> const& mat) const;
		/*!Multiplies two matrices (m1 *= m2 : m1 = m1*m2)*/
		Matrix<Type>& operator*=(Matrix<Type> const& mat);
		/*!Calls operator*=(Matrix<Type> const& mat)*/
		Matrix<Type> operator*(Matrix<Type> const& mat) const;

		/*!Multiplies a matrix by a scalar Type*/
		Matrix<Type>& operator*=(Type const& d);
		/*!Calls operator*=(Type const& d)*/
		Matrix<Type> operator*(Type const& d) const;
		/*!Devides a matrix by a scalar*/
		Matrix<Type>& operator/=(Type const& d);
		/*!Substracts a matrix by a scalar*/
		Matrix<Type>& operator-=(Type const& d);
		/*!Add a matrix by a scalar*/
		Matrix<Type>& operator+=(Type const& d);

		/*!Multiplies a matrix by a vector (m1 *= v2 : m1 = m1*v2)*/
		Vector<Type> operator*(Vector<Type> const& vec) const;

		/*!Set the matrix to 0*/
		void set();
		/*!Set the matrix to val*/
		void set(Type const& val);
		/*!Set a N_row x N_col matrix  */
		void set(unsigned int N_row, unsigned int N_col);
		/*!Set a N_row x N_col matrix to val */
		void set(unsigned int N_row, unsigned int N_col, Type val);

		/*!Sets the entries to zero if they are close to 0*/
		Matrix<Type> chop(double precision = 1e-10) const;

		/*!Print the matrix for mathematica*/
		void print_mathematica() const;

		/*!Returns the transpose of any matrix*/
		Matrix<Type> transpose() const;
		/*!Returns the conjugate transpose of complex matrix*/
		Matrix<Type> conjugate_transpose() const;
		/*!Returns the diagonal elements in an vector*/
		Vector<Type> diag() const;
		/*!Returns the trace*/
		Type trace() const;

		/*!Returns the size of the matrix*/
		unsigned int const& size() const { return size_; }
		/*!Returns the number of rows of the matrix*/
		unsigned int const& row() const { return row_; }
		/*!Returns the number of columns of the matrix*/
		unsigned int const& col() const { return col_; }
		/*!Returns the pointer to the matrix*/
		Type* ptr() const { return mat_; }

		/*!Returns the maximal value of mat_*/
		Type max() const;

		std::string header_def() const { return "Matrix("+RST::math(my::tostring(row_)+"\\times"+my::tostring(col_))+")"; }

	protected:
		unsigned int row_  = 0; //!< number of rows
		unsigned int col_  = 0; //!< number of columns
		unsigned int size_ = 0; //!< size of the array
		Type* mat_ = NULL;		//!< pointer to a static array

		/*!Copy-And-Swap Idiom*/
		void swap_to_assign(Matrix<Type>& m1,Matrix<Type>& m2);
};

/*constructors and destructor*/
/*{*/
template<typename Type>
Matrix<Type>::Matrix(unsigned int N_row, unsigned int N_col):
	row_(N_row),
	col_(N_col),
	size_(N_col*N_row),
	mat_(size_?new Type[size_]:NULL)
{}

template<typename Type>
Matrix<Type>::Matrix(unsigned int N_row, unsigned int N_col, Type val):
	row_(N_row),
	col_(N_col),
	size_(N_col*N_row),
	mat_(size_?new Type[size_]:NULL)
{
	for(unsigned int i(0);i<size_;i++){ mat_[i] = val; }
}

template<typename Type>
Matrix<Type>::Matrix(Matrix<Type> const& mat):
	row_(mat.row_),
	col_(mat.col_),
	size_(mat.size_),
	mat_(size_?new Type[size_]:NULL)
{
	for(unsigned int i(0);i<size_;i++){ mat_[i] = mat.mat_[i]; }
}

template<typename Type>
Matrix<Type>::Matrix(Matrix<Type>&& mat):
	row_(mat.row_),
	col_(mat.col_),
	size_(mat.size_),
	mat_(mat.mat_)
{
	mat.mat_ = NULL;
	mat.row_ = mat.col_ = mat.size_ = 0;
}

template<typename Type>
Matrix<Type>::Matrix(IOFiles& r):
	row_(r.read<unsigned int>()),
	col_(r.read<unsigned int>()),
	size_(row_*col_),
	mat_(size_?new Type[size_]:NULL)
{
	r.read(mat_,size_,sizeof(Type));
}

template<typename Type>
Matrix<Type>::~Matrix(){
	if(mat_){
		delete[]  mat_;
		mat_ = NULL;
	}
}

template<typename Type>
void Matrix<Type>::swap_to_assign(Matrix<Type>& m1,Matrix<Type>& m2){
	std::swap(m1.mat_,m2.mat_);
	std::swap(m1.row_,m2.row_);
	std::swap(m1.col_,m2.col_);
	std::swap(m1.size_,m2.size_);
}
/*}*/

/*i/o methods*/
/*{*/
template<typename Type>
std::ostream& operator<<(std::ostream& flux, Matrix<Type> const& m){
	for(unsigned int i(0);i<m.row();i++){
		for(unsigned int j(0);j<m.col();j++){ flux<<m(i,j)<<" "; }
		if(i+1 != m.row()){ flux<<std::endl; }
	}
	return flux;
}

template<typename Type>
std::istream& operator>>(std::istream& flux, Matrix<Type>& m){
	unsigned int row(m.row());
	unsigned int col(m.col());
	for(unsigned int i(0);i<row;i++){
		for(unsigned int j(0);j<col;j++){
			flux>>m.ptr()[i+j*row];
		}
	}
	return flux;
}

template<typename Type>
IOFiles& operator<<(IOFiles& w, Matrix<Type> const& m){
	if(w.is_binary()){
		w<<m.row()<<m.col();
		w.write(m.ptr(),m.size(),sizeof(Type));
	} else { w.stream()<<m; }
	return w;
}

template<typename Type>
IOFiles& operator>>(IOFiles& r, Matrix<Type>& m){
	if(r.is_binary()){ m = std::move(Matrix<Type>(r)); }
	else { r.stream()>>m; }
	return r;
}
/*}*/

/*arithmetic operators*/
/*{*/
/*{Matrix.operator(Matrix)*/
template<typename Type>
Matrix<Type>& Matrix<Type>::operator=(Matrix<Type> mat){
	swap_to_assign(*this,mat);
	return (*this);
}

template<typename Type>
Matrix<Type>& Matrix<Type>::operator+=(Matrix<Type> const& mat){
	assert(row_ == mat.row_ && col_ == mat.col_);
	for(unsigned int i(0);i<size_;i++){ mat_[i] += mat.mat_[i]; }
	return (*this);
}

template<typename Type>
Matrix<Type> Matrix<Type>::operator+(Matrix<Type> const& mat) const {
	Matrix<Type> matout((*this));
	matout += mat;
	return matout;
}

template<typename Type>
Matrix<Type>& Matrix<Type>::operator-=(Matrix<Type> const& mat){
	assert(row_ == mat.row_ && col_ == mat.col_);
	for(unsigned int i(0);i<size_;i++){ mat_[i] -= mat.mat_[i]; }
	return (*this);
}

template<typename Type>
Matrix<Type> Matrix<Type>::operator-(Matrix<Type> const& mat) const {
	Matrix<Type> matout((*this));
	matout -= mat;
	return matout;
}

template<typename Type>
Matrix<Type>& Matrix<Type>::operator*=(Matrix<Type> const& mat){
	assert(col_==mat.row_);
	Matrix<Type> tmp(*this);
	if(col_!=row_ || mat.col_!=mat.row_ ){
		if(mat_){ delete[] mat_; }
		mat_ = new Type[row_*mat.col_];
		size_ = row_*mat.col_;
		col_ = mat.col_;
	}
	for(unsigned int i(0);i<row_;i++){
		for(unsigned int j(0);j<col_;j++){
			mat_[i+j*row_] = 0.0;
			for(unsigned int k(0);k<tmp.col_;k++){
				mat_[i+j*row_] += tmp.mat_[i+k*tmp.row_] * mat.mat_[k+j*mat.row_];
			}
		}
	}
	return (*this);
}

template<typename Type>
Matrix<Type> Matrix<Type>::operator*(Matrix<Type> const& mat) const {
	assert(col_==mat.row_);
	Matrix<Type> matout(row_,mat.col_);
	for(unsigned int i(0);i<matout.row_;i++){
		for(unsigned int j(0);j<matout.col_;j++){
			matout.mat_[i+j*matout.row_] = 0.0;
			for(unsigned int k(0);k<mat.row_;k++){
				matout.mat_[i+j*matout.row_] += mat_[i+k*row_] * mat.mat_[k+j*mat.row_];
			}
		}
	}
	return matout;
}

template<typename Type>
Matrix<Type> operator*(Type const& d, Matrix<Type> const& mat) {
	return mat*d;
}
/*}*/

/*{Matrix.operator(Type)*/
template<typename Type>
Matrix<Type>& Matrix<Type>::operator*=(Type const& d){
	for(unsigned int i(0);i<size_;i++){
		mat_[i] *= d;
	}
	return (*this);
}

template<typename Type>
Matrix<Type> Matrix<Type>::operator*(Type const& d) const {
	Matrix<Type> tmp(*this);
	tmp *= d;
	return tmp;
}

template<typename Type>
Matrix<Type>& Matrix<Type>::operator-=(Type const& d){
	for(unsigned int i(0);i<size_;i++){ mat_[i] -= d; }
	return (*this);
}

template<typename Type>
Matrix<Type>& Matrix<Type>::operator+=(Type const& d){
	for(unsigned int i(0);i<size_;i++){ mat_[i] += d; }
	return (*this);
}

template<typename Type>
Matrix<Type>& Matrix<Type>::operator/=(Type const& d){
	for(unsigned int i(0);i<size_;i++){ mat_[i] /= d; }
	return (*this);
}
/*}*/

/*{Matrix.operator(Vector)*/
template<typename Type>
Vector<Type> Matrix<Type>::operator*(Vector<Type> const& vec) const {
	assert(vec.size() == col_);
	Vector<Type> tmp(row_);
	for(unsigned int i(0);i<row_;i++){
		tmp(i) = 0.0;
		for(unsigned int j(0);j<col_;j++){
			tmp(i) += mat_[i+j*row_] * vec(j);
		}
	}
	return tmp;
}
/*}*/
/*}*/

/*methods that modify the class*/
/*{*/
template<>
inline Matrix<double> Matrix<double>::chop(double precision) const {
	Matrix<double> tmp(*this);
	for(unsigned int i(0);i<tmp.size_;i++){
		if(std::abs(tmp.mat_[i]) < precision ){ tmp.mat_[i]=0.0; }
	}
	return tmp;
}

template<>
inline Matrix<std::complex<double> > Matrix<std::complex<double> >::chop(double precision) const {
	Matrix<std::complex<double> > tmp(*this);
	for(unsigned int i(0);i<tmp.size_;i++){
		if(std::abs(tmp.mat_[i].imag()) < precision ){ tmp.mat_[i].imag(0.0); }
		if(std::abs(tmp.mat_[i].real()) < precision ){ tmp.mat_[i].real(0.0); }
	}
	return tmp;
}

template<typename Type>
void Matrix<Type>::set(){
	if(mat_){ delete[] mat_; }
	mat_ = NULL;
	row_ = 0;
	col_ = 0;
	size_ = 0;
}

template<typename Type>
void Matrix<Type>::set(Type const& val){
	for(unsigned int i(0); i<size_; i++){
		mat_[i] = val;
	}
}

template<typename Type>
void Matrix<Type>::set(unsigned int row, unsigned int col){
	if(col_ != col || row_ != row){
		if(mat_){ delete[] mat_; }
		mat_ = new Type[row*col];
		row_ = row;
		col_ = col;
		size_ = row*col;
	}
}

template<typename Type>
void Matrix<Type>::set(unsigned int row, unsigned int col, Type val){
	this->set(row,col);
	this->set(val);
}
/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename Type>
Matrix<Type> Matrix<Type>::transpose() const {
	Matrix<Type> tmp(col_,row_);
	for(unsigned int i(0);i<col_;i++){
		for(unsigned int j(0);j<row_;j++){
			tmp.mat_[i+j*tmp.row_] = mat_[j+i*row_];
		}
	}
	return tmp;
}

template<>
inline Matrix<std::complex<double> > Matrix<std::complex<double> >::conjugate_transpose() const {
	Matrix<std::complex<double> > tmp(col_,row_);
	for(unsigned int i(0);i<col_;i++){
		for(unsigned int j(0);j<row_;j++){
			tmp.mat_[i+j*tmp.row_] = std::conj(mat_[j+i*row_]);
		}
	}
	return tmp;
}

template<>
inline Matrix<double> Matrix<double>::conjugate_transpose() const {
	std::cerr<<__PRETTY_FUNCTION__<<" : not defined for real matrices, replaced by transpose()"<<std::endl;
	return (*this).transpose();
}

template<typename Type>
Vector<Type> Matrix<Type>::diag() const {
	unsigned int N(0);
	if(row_<col_){
		N=col_;
		std::cerr<<__PRETTY_FUNCTION__<<" : to check"<<std::endl;
	} else { N=row_; }
	Vector<Type> v(N);
	for(unsigned int i(0);i<N;i++){
		v(i) = mat_[i*(row_+1)];
	}
	return v;
}

template<>
inline void Matrix<double>::print_mathematica() const {
	std::cout<<"{{";
	for(unsigned int i(0);i<row_;i++){
		for(unsigned int j(0);j<col_;j++){
			std::cout<<mat_[i+j*row_];
			if(j+1==col_){
				if(i+1==row_){ std::cout<<"}}"<<std::endl; }
				else{ std::cout<<"},"<<std::endl<<"{"; }
			} else{ std::cout<<","; }
		}
	}
}

template<>
inline void Matrix<std::complex<double> >::print_mathematica() const {
	std::complex<double> a;
	std::cout<<"{";
	for(unsigned int i(0);i<row_;i++){
		std::cout<<"{ ";
		for(unsigned int j(0);j<col_;j++){
			a=mat_[i+j*row_];
			std::cout<<a.real();
			if(a.imag()>0){ std::cout<<"+"; }
			else{ std::cout<<"-"; }
			std::cout<<a.imag()<<" I ";
			if(j+1==col_){
				if(i+1==row_){ std::cout<<"}}"<<std::endl; }
				else {std::cout<<"},"<<std::endl<<"{"; }
			} else { std::cout<<","; }
		}
	}
}

template<typename Type>
Type Matrix<Type>::trace() const {
	unsigned int k(std::min(row_,col_));
	Type t(0.0);
	std::cerr<<__PRETTY_FUNCTION__<<" : need to be checked"<<std::endl;
	for(unsigned int i(0);i<k;i++){
		t += mat_[i*(row_+1)];
	}
	return t;
}

template<typename Type>
Type Matrix<Type>::max() const {
	Type m(mat_[0]);
	for(unsigned int i(1);i<size_;i++){
		if(mat_[i]>m){ m=mat_[i]; }
	}
	return m;
}
/*}*/
#endif
