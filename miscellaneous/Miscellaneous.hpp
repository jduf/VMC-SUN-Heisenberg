#ifndef DEF_MISCELLANEOUS
#define DEF_MISCELLANEOUS

#include <string>
#include <iostream>
#include <sstream>
template<typename Type>
std::string tostring(Type const& t){
	std::ostringstream s;
	s<<t;
	return s.str();
}

#include <complex>
#include <cmath> 
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

/*double chop(T)*/
/*{*/
inline double chop(double const& x, double precision = 1e-10){ return (std::abs(x)<precision?0.0:x); }

inline std::complex<double> chop(std::complex<double> x, double precision = 1e-10){
	if(std::abs(x.imag()) < precision ){x.imag(0.0);}
	if(std::abs(x.real()) < precision ){x.real(0.0);}
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

#ifdef DEF_VECTOR
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
#endif
/*}*/
#endif
