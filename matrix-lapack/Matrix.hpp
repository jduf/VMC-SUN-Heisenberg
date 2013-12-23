#ifndef DEF_MATRIX
#define DEF_MATRIX

#include <iostream>
#include <cassert>
#include <cmath> //allow abs(double) and abs(complex) 
#include <complex>

//#include "Vector.hpp"

template<typename Type>
class Vector;

/*!Class that implement a static array as a Matrix
 *
 * - can be saved with Write.hpp 
 * - can be loaded with Read.hpp 
*/
template<typename Type>
class Matrix{
	public:
		/*!Default constructor that initializes *m to NULL and N to 0*/
		Matrix();
		/*!Initializes a static array of Type of size N_row*N_col*/
		Matrix(unsigned int N_row, unsigned int N_col);
		/*!Initializes a static array of Type of size N_row*N_col to a value val*/
		Matrix(unsigned int N_row, unsigned int N_col, Type val);
		/*!Deep copy*/
		Matrix(Matrix<Type> const& mat);
		/*!Replace a pointer with an instance*/
		Matrix(Matrix<Type> *map);
		/*!Delete the static array*/
		~Matrix();

		/*!Accesses the (i,j)th entry of the matrix*/
		Type const& operator()(unsigned int const& i, unsigned int const& j)
			const { assert(i<row_ && j<col_); return m_[i+j*row_]; };
		/*!Sets the (i,j)th entry of the matrix*/
		Type& operator()(unsigned int const& i, unsigned int const& j) {
			assert(i<row_ && j<col_); return m_[i+j*row_]; };

		/*!Deep copy assignment*/
		Matrix<Type>& operator=(Matrix<Type> const& mat); 
		/*!Additions this matrice with another*/
		Matrix<Type>& operator+=(Matrix<Type> const& mat);
		Matrix<Type> operator+(Matrix<Type> const& mat) const;
		/*!Substracts this matrice from another (m1 -= m2 : m1 = m1-m2)*/
		Matrix<Type>& operator-=(Matrix<Type> const& mat);
		Matrix<Type> operator-(Matrix<Type> const& mat) const;
		/*!Multiplies two matrices (m1 *= m2 : m1 = m1*m2)*/
		Matrix<Type>& operator*=(Matrix<Type> const& mat);
		Matrix<Type> operator*(Matrix<Type> const& mat) const;

		/*!Multiplies a matrix by a scalar Type*/
		Matrix<Type>& operator*=(Type const& d);
		Matrix<Type> operator*(Type const& d) const;

		Matrix<Type>& operator-=(Type const& d);
		Matrix<Type>& operator/=(Type const& d);

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
		/*!Returns the conjugate transpose of complex matrix (may give an error) */
		Matrix<Type> trans_conj() const;
		/*!Returns the diagonal elements in an vector*/
		Vector<Type> diag() const;
		/*!Returns the trace*/
		Type trace() const;

		/*!Returns the pointer to the matrix*/
		Type* ptr() const { return m_; }
		/*!Returns the size of the matrix*/
		unsigned int size() const { return size_; }
		/*!Returns the number of rows of the matrix*/
		unsigned int row() const { return row_; }
		/*!Returns the number of columns of the matrix*/
		unsigned int col() const { return col_; }

		/*return the maximum value of the matrix*/
		Type max() const;

	protected:
		Type *m_; //!< pointer to a static array
		unsigned int row_; //!< number of rows
		unsigned int col_; //!< number of columns
		unsigned int size_; //!< size of the array

		void set_null_pointer(){m_=NULL;}
};

template<typename Type>
std::ostream& operator<<(std::ostream& flux, Matrix<Type> const& mat);
template<typename Type>
std::istream& operator>>(std::istream& flux, Matrix<Type>& mat);

/*constructors and destructor*/
/*{*/
template<typename Type>
Matrix<Type>::Matrix():
	m_(NULL),
	row_(0),
	col_(0),
	size_(0)
{ }

template<typename Type>
Matrix<Type>::Matrix(unsigned int N_row, unsigned int N_col):
	m_(new Type[N_row*N_col]),
	row_(N_row),
	col_(N_col),
	size_(N_col*N_row)
{ } 

template<typename Type>
Matrix<Type>::Matrix(unsigned int N_row, unsigned int N_col, Type val):
	m_(new Type[N_row*N_col]),
	row_(N_row),
	col_(N_col),
	size_(N_col*N_row)
{ 
	set(val);
}

template<typename Type>
Matrix<Type>::Matrix(Matrix<Type> const& mat):
	m_(new Type[mat.size_]),
	row_(mat.row_),
	col_(mat.col_),
	size_(mat.size_)
{
	for(unsigned int i(0);i<size_;i++){
		m_[i] = mat.m_[i];
	}
}

template<typename Type>
Matrix<Type>::Matrix(Matrix<Type> *mat):
	m_(mat->m_),
	row_(mat->row_),
	col_(mat->col_),
	size_(mat->size_)
{ 
	mat->set_null_pointer();
}

template<typename Type>
Matrix<Type>::~Matrix(){
	if(m_){
		delete[]  m_;
		m_ = NULL;
	}
}
/*}*/

/*operators*/
/*{*/
template<typename Type>
Matrix<Type>& Matrix<Type>::operator=(Matrix<Type> const& mat){
	if(!mat.m_){ 
		set();
	} else {
		if(col_ != mat.col_ ||  row_ != mat.row_){
			if(m_){ delete[] m_;}
			m_ = new Type[mat.size_];
			size_ = mat.size_;
			row_ = mat.row_;
			col_ = mat.col_;
		}
		for(unsigned int i(0); i<size_; i++){
			m_[i] = mat.m_[i];
		}
	}
	return (*this);
}

template<typename Type>
Matrix<Type>& Matrix<Type>::operator+=(Matrix<Type> const& mat){
	assert(row_ == mat.row_ && col_ == mat.col_);
	for(unsigned int i(0);i<size_;i++){
		m_[i] += mat.m_[i];
	}
	return (*this);
}

template<typename Type>
Matrix<Type> Matrix<Type>::operator+(Matrix<Type> const& mat) const{
	Matrix<Type> matout((*this));
	matout += mat;
	return matout;
}

template<typename Type>
Matrix<Type>& Matrix<Type>::operator-=(Matrix<Type> const& mat){
	assert(row_ == mat.row_ && col_ == mat.col_);
	for(unsigned int i(0);i<size_;i++){
		m_[i] -= mat.m_[i];
	}
	return (*this);
}

template<typename Type>
Matrix<Type> Matrix<Type>::operator-(Matrix<Type> const& mat) const{
	Matrix<Type> matout((*this));
	matout -= mat;
	return matout;
}

template<typename Type>
Matrix<Type>& Matrix<Type>::operator*=(Matrix<Type> const& mat){
	assert(col_==mat.row_);
	Matrix<Type> tmp(*this);
	if(col_!=row_ || mat.col_!=mat.row_ ){
		if(m_){ delete[] m_;}
		m_ = new Type[row_*mat.col_];
		size_ = row_*mat.col_;
		col_ = mat.col_;
	}
	for(unsigned int i(0);i<row_;i++){
		for(unsigned int j(0);j<col_;j++){
			m_[i+j*row_] = 0.0;
			for(unsigned int k(0);k<tmp.col_;k++){
				m_[i+j*row_] += tmp.m_[i+k*tmp.row_] * mat.m_[k+j*mat.row_];
			}
		}
	}
	return (*this);
}

template<typename Type>
Matrix<Type> Matrix<Type>::operator*(Matrix<Type> const& mat) const{
	assert(col_==mat.row_);
	Matrix<Type> matout(row_,mat.col_);
	for(unsigned int i(0);i<matout.row_;i++){
		for(unsigned int j(0);j<matout.col_;j++){
			matout.m_[i+j*matout.row_] = 0.0;
			for(unsigned int k(0);k<mat.row_;k++){
				matout.m_[i+j*matout.row_] += m_[i+k*row_] * mat.m_[k+j*mat.row_];
			}
		}
	}
	return matout;
}

template<typename Type>
std::ostream& operator<<(std::ostream& flux, Matrix<Type> const& mat){
	for(unsigned int i(0);i<mat.row();i++){
		for(unsigned int j(0);j<mat.col();j++){
			flux<<mat(i,j)<<" ";
		}
		if(i+1 != mat.row()){ flux<<std::endl; }
	}
	return flux;
}

template<typename Type>
std::istream& operator>>(std::istream& flux, Matrix<Type>& mat){
	for(unsigned int i(0);i<mat.row();i++){
		for(unsigned int j(0);j<mat.col();j++){
			flux>>mat(i,j);
		}
	}
	return flux;
}

template<typename Type>
Matrix<Type>& Matrix<Type>::operator*=(Type const& d){
	for(unsigned int i(0);i<size_;i++){
		m_[i] *= d; 
	}
	return (*this);
}

template<typename Type>
Matrix<Type> Matrix<Type>::operator*(Type const& d) const{
	Matrix<Type> tmp(*this);
	tmp *= d;
	return tmp;
}

template<typename Type>
Matrix<Type> operator*(Type const& d, Matrix<Type> const& mat) {
	return mat*d;
}

template<typename Type>
Matrix<Type>& Matrix<Type>::operator-=(Type const& d){
	for(unsigned int i(0);i<size_;i++){
		m_[i] -= d; 
	}
	return (*this);
}

template<typename Type>
Matrix<Type>& Matrix<Type>::operator/=(Type const& d){
	for(unsigned int i(0);i<size_;i++){
		m_[i] /= d; 
	}
	return (*this);
}
/*}*/

/*methods that modify the class*/
/*{*/
template<>
inline Matrix<double> Matrix<double>::chop(double precision) const {
	Matrix<double> tmp(*this);
	for(unsigned int i(0);i<tmp.size_;i++){
		if(std::abs(tmp.m_[i]) < precision ){tmp.m_[i]=0.0;}
	}
	return tmp;
}

template<>
inline Matrix<std::complex<double> > Matrix<std::complex<double> >::chop(double precision) const{
	Matrix<std::complex<double> > tmp(*this);
	for(unsigned int i(0);i<tmp.size_;i++){
		if(std::abs(tmp.m_[i].imag()) < precision ){tmp.m_[i].imag(0.0);}
		if(std::abs(tmp.m_[i].real()) < precision ){tmp.m_[i].real(0.0);}
	}
	return tmp;
}

template<typename Type>
void Matrix<Type>::set(){
	if(m_){ delete[] m_; }
	m_ = NULL;
	row_ = 0;
	col_ = 0;
	size_ = 0;
}

template<typename Type>
void Matrix<Type>::set(Type const& val){
	for(unsigned int i(0); i<size_; i++){
		m_[i] = val;
	}
}

template<typename Type>
void Matrix<Type>::set(unsigned int row, unsigned int col){
	if(col_ != col || row_ != row){ 
		if(m_){ delete[] m_; }
		m_ = new Type[row*col];
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

//template<typename Type>
//void Matrix<Type>::apply_index(Matrix<unsigned int> const& index){
//Matrix<Type> tmp(*this);
//for(unsigned int i(0);i<row_;i++){
//(*this)(i) = tmp(index(i));
//}
//}

/*}*/

/*methods that return something related to the class*/
/*{*/
template<typename Type>
Matrix<Type> Matrix<Type>::transpose() const{
	Matrix<Type> tmp(col_,row_);
	for(unsigned int i(0);i<col_;i++){
		for(unsigned int j(0);j<row_;j++){
			tmp.m_[i+j*tmp.row_] = m_[j+i*row_];
		}
	}
	return tmp;
}

template<>
inline Matrix<std::complex<double> > Matrix<std::complex<double> >::trans_conj() const{
	Matrix<std::complex<double> > tmp(col_,row_);
	for(unsigned int i(0);i<col_;i++){
		for(unsigned int j(0);j<row_;j++){
			tmp.m_[i+j*tmp.row_] = std::conj(m_[j+i*row_]);
		}
	}
	return tmp;
}

template<>
inline Matrix<double> Matrix<double>::trans_conj() const {
	std::cerr<<"Matrix : conj_trans : not defined for real matrices"<<std::endl;
	std::cerr<<"Matrix : conj_trans : replace by transpose()"<<std::endl;
	return (*this).transpose();
}

template<typename Type>
Vector<Type> Matrix<Type>::diag() const{
	unsigned int N(0);
	if(row_ < col_){N=col_; std::cerr<<"Matrix : diag : to check"<<std::endl; }
	else{N=row_;}
	Vector<Type> v(N);
	for(unsigned int i(0);i<N;i++){
		v(i) = m_[i*(row_+1)];
	}
	return v;
}

template<>
inline void Matrix<double>::print_mathematica() const {
	std::cout<<"{{";
	for(unsigned int i(0);i<row_;i++){
		for(unsigned int j(0);j<col_;j++){
			std::cout<<m_[i+j*row_];
			if(j+1==col_){
				if(i+1==row_){std::cout<<"}}"<<std::endl;}
				else{std::cout<<"},"<<std::endl<<"{";}
			} else{
				std::cout<<",";
			}
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
			a=m_[i+j*row_];
			std::cout<<a.real();
			if(a.imag()>0){std::cout<<"+";}
			else{std::cout<<"-";}
			std::cout<<a.imag()<<" I ";
			if(j+1==col_){
				if(i+1==row_){std::cout<<"}}"<<std::endl;}
				else{std::cout<<"},"<<std::endl<<"{";}
			} else{
				std::cout<<",";
			}
		}
	}
}

template<typename Type>
Type Matrix<Type>::trace() const {
	unsigned int k(std::min(row_,col_));
	Type t(0.0);
	std::cerr<<"Matrix : trace : need to be checked"<<std::endl;
	for(unsigned int i(0);i<k;i++){
		t += m_[i*(row_+1)];
	}
	return t;
}

template<typename Type>
Type Matrix<Type>::max() const {
	Type m(m_[0]);
	for(unsigned int i(1);i<size_;i++){
		if(m_[i]>m){m=m_[i];}
	}
	return m;
}
/*}*/

/*double real(T)*/
/*{*/
inline double real(double const& x){ return x; }

inline double real(std::complex<double> const& x){ return std::real(x); }
/*}*/

/*double norm_squared(T)*/
/*{*/
inline double norm_squared(double x){
	return x*x;
}

inline double norm_squared(std::complex<double> x){
	return std::norm(x);
}
/*}*/
#endif
