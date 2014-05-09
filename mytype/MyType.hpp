#ifndef DEF_MYTYPE
#define DEF_MYTYPE

#include <complex>
#include <iostream>
#include <cmath> //allow abs(double) and abs(complex) 
#include <complex>
#include <cassert>

#ifdef DEF_IOFILES
class IOFiles;
class RST;
#endif

template<typename Type>
class MyType{
	public:
		virtual void write_in_stream(std::ostream& flux) const = 0;
		virtual void read_from_stream(std::istream& flux) = 0;
#ifdef DEF_IOFILES
		virtual void write_in_file(IOFiles& w) const = 0;
		virtual void read_from_file(IOFiles& r) = 0;
		virtual void header_rst(std::string const& s, RST& rst) const { rst.title(s,"=");}
#endif
};

template<typename Type>
std::ostream& operator<<(std::ostream& flux, MyType<Type> const& t){
	t.write_in_stream(flux);
	return flux;
}

template<typename Type>
std::istream& operator>>(std::istream& flux, MyType<Type>& t){
	t.read_from_stream(flux);
	return flux;
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
inline double norm_squared(double x){
	return x*x;
}

inline double norm_squared(std::complex<double> x){
	return std::norm(x);
}
/*}*/
#endif
