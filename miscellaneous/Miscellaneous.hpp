#ifndef DEF_MISCELLANEOUS
#define DEF_MISCELLANEOUS

#include <string>
#include <iostream>
#include <sstream>
#include <complex>
#include <cmath>
#include <vector>
#include <algorithm>

namespace my{
	template<typename Type>
		std::string tostring(Type const& t){
			std::ostringstream s;
			s<<t;
			return s.str();
		}

	template<typename Type>
		bool string2type(std::string const& s, Type& out){
			std::stringstream ss(s);
			if( ss>>out ){ return true; }
			else {
				std::cerr<<__PRETTY_FUNCTION__<<" incoherent types"<<std::endl;
				return false;
			}
		}

	template<typename Type>
		Type string2type(std::string const& s){
			Type t;
			my::string2type(s,t);
			return t;
		}

	inline double get_double(std::string const& msg){
		std::string token;
		double in;
		do{
			std::cout<<msg<<" [double] ";
			std::getline(std::cin,token);
		} while (!(token.find(' ') == std::string::npos && my::string2type<double>(token,in)));
		return in;
	}

	inline bool get_yn(std::string const& msg){
		std::string token;
		do{
			std::cout<<msg<<" [y/n] ";
			std::getline(std::cin,token);
			if(token == "y"){ return true; }
			if(token == "n"){ return false; }
		} while(1);
	}

	inline std::vector<std::string> &string_split(const std::string &s, char delim, std::vector<std::string> &elems){
		std::stringstream ss(s);
		std::string item;
		while (std::getline(ss, item, delim)) { elems.push_back(item); }
		return elems;
	}

	inline std::vector<std::string> string_split(const std::string &s, char delim){
		std::vector<std::string> elems;
		string_split(s, delim, elems);
		if(elems.back()==""){ elems.pop_back(); }
		return elems;
	}

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
	inline double norm_squared(double const& x){ return x*x; }

	inline double norm_squared(std::complex<double> const& x){ return std::norm(x); }
	/*}*/

	/*double chop(T)*/
	/*{*/
	inline double chop(double const& x, double precision = 1e-10){ return (std::abs(x)<precision?0.0:x); }

	inline std::complex<double> chop(std::complex<double> x, double precision = 1e-10){
		if(std::abs(x.imag()) < precision ){ x.imag(0.0); }
		if(std::abs(x.real()) < precision ){ x.real(0.0); }
		return x;
	}
	/*}*/

	/*bool are_equal(T,T)*/
	/*{*/
	inline bool are_equal(double x, double y, double abs_tol=1e-14, double rel_tol=1e-14){
		double diff(std::abs(x-y));
		x = std::abs(x);
		y = std::abs(y);
		x = (x>y)?x:y;
		return (diff<rel_tol*x) || (diff<abs_tol);
	}

	inline bool are_equal(std::complex<double> const& x, std::complex<double> const& y, double abs_tol=1e-14, double rel_tol=1e-14){
		if(!are_equal(std::abs(x),std::abs(y),abs_tol,rel_tol)){ return false; }
		if(!are_equal(x.real(),y.real(),abs_tol,rel_tol)){ return false; }
		if(!are_equal(x.imag(),y.imag(),abs_tol,rel_tol)){ return false; }
		return true;
	}
	/*}*/

	inline unsigned long long gcd(unsigned long long x, unsigned long long y){
		unsigned long long t;
		while(y != 0){
			t = x % y;
			x = y;
			y = t;
		}
		return x;
	}

	inline unsigned long long nCk(unsigned long long n, unsigned long long k){
		if(k > n){ std::cerr<<__PRETTY_FUNCTION__<<" : k > n"<<std::endl; }
		unsigned long long r(1);
		unsigned long long g;
		unsigned long long t;
		for (unsigned long long d(1); d<=k;++d,--n){
			g = gcd(r,d);
			r /= g;
			t = n / (d / g);
			if(r > std::numeric_limits<unsigned long long>::max() / t){
				std::cerr<<__PRETTY_FUNCTION__<<" : other bug"<<std::endl;
			}
			r *= t;
		}
		return r;
	}
}

namespace BLAS{
	/*{BLAS level 1*/
	extern "C" double ddot_(unsigned int const& N, double const* const dx, unsigned int const& ix, double const* const dy, unsigned int const& iy);
	inline double dot(
			unsigned int const& N,
			double const* const a,
			bool const& ar, //true : multiply a row of a
			unsigned int const& arow,// 1 for Vector 
			unsigned int aidx, // 0 for Vector
			double const* const b, 
			bool const& br, //true : multiply a row of b
			unsigned int const& brow,// 1 for Vector 
			unsigned int bidx// 0 for Vector
			)
	{
		if(ar){//multiply a row of a
			if(br){ return ddot_(N,/**/a+aidx,     arow,/**/b+bidx     ,brow); }
			else  { return ddot_(N,/**/a+aidx,     arow,/**/b+bidx*brow,1); } 
		} else {                                         
			std::cout<<__PRETTY_FUNCTION__<<" : need to be checked"<<std::endl; 
			if(br){ return ddot_(N,/**/a+aidx*arow,1,   /**/b+bidx     ,brow); }
			else  { return ddot_(N,/**/a+aidx*arow,1,   /**/b+bidx*brow,1); }
		}
	}

	extern "C" std::complex<double> zdotu_(unsigned int const& N, std::complex<double> const* const dx, unsigned int const& ix, std::complex<double> const* const dy, unsigned int const& iy);
	inline std::complex<double> dot(
			unsigned int const& N,
			std::complex<double> const* const a,
			bool const& ar, //true : multiply a row of a
			unsigned int const& arow,// 1 for Vector 
			unsigned int aidx, // 0 for Vector
			std::complex<double> const* const b, 
			bool const& br,  //true : multiply a row of b
			unsigned int const& brow,// 0 for Vector 
			unsigned int bidx// 1 for Vector
			)
	{
		if(ar){
			if(br){ return zdotu_(N,/**/a+aidx,     arow,/**/b+bidx     ,brow); }
			else  { return zdotu_(N,/**/a+aidx,     arow,/**/b+bidx*brow,1); } 
		} else {
			std::cout<<__PRETTY_FUNCTION__<<" : need to be checked"<<std::endl; 
			if(br){ return zdotu_(N,/**/a+aidx*arow,1,   /**/b+bidx     ,brow); }
			else  { return zdotu_(N,/**/a+aidx*arow,1,   /**/b+bidx*brow,1); }
		}
	}
	/*}*/

	/*{BLAS level 2*/
	extern "C" void dgemv_(char const& trans, unsigned int const& N, unsigned int const& M, double const& alpha, double const* const a, unsigned int const& lda, double const* const x, unsigned int const& incx, double const& beta, double const* y, unsigned int const& incy);
	inline void gemv(
			char const& trans,
			unsigned int const& N,
			unsigned int const& M,
			double const* const a,
			double const* const x,
			unsigned int const& incx,
			double const* y
			)
	{
		dgemv_(trans,N,M,1.0,a,M,x,incx,0.0,y,1);
	}

	extern "C" void zgemv_(char const& trans, unsigned int const& N, unsigned int const& M, std::complex<double> const& alpha, std::complex<double> const* const a, unsigned int const& lda, std::complex<double> const* const x, unsigned int const& incx, std::complex<double> const& beta, std::complex<double> const* y, unsigned int const& incy);
	inline void gemv(
			char const& trans,
			unsigned int const& N,
			unsigned int const& M,
			std::complex<double> const* const a,
			std::complex<double> const* const x,
			unsigned int const& incx,
			std::complex<double> const* y
			)
	{
		zgemv_(trans,N,M,1.0,a,M,x,incx,0.0,y,1);
	}
	/*}*/
}
#endif
