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

	inline std::vector<std::string> &string_split(std::string const& s, char const& delim, std::vector<std::string>& elems){
		std::stringstream ss(s);
		std::string item;
		while (std::getline(ss, item, delim)) { elems.push_back(item); }
		return elems;
	}

	inline std::vector<std::string> string_split(std::string const& s, char const& delim){
		std::vector<std::string> elems;
		string_split(s, delim, elems);
		if(elems.back()==""){ elems.pop_back(); }
		return elems;
	}

	inline void ensure_trailing_slash(std::string& s){ if(s.back() != '/'){ s += '/'; } }

	/*int to alphabet*/
	inline char int_to_alphabet(unsigned int const& i, bool const& upper_case){
		if(i<26){
			if(upper_case){ return "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[i]; }
			else { return "abcdefghijklmnopqrstuvwxyz"[i]; }
		} else { return '?'; }
	}

	/*round*/
	/*{*/
	inline double round_nearest(double const& d, unsigned int const decimal_place){
		return roundf(d*decimal_place)/decimal_place;
	}

	inline double round_down(double const& d, unsigned int const decimal_place){
		return floorf(d*decimal_place)/decimal_place;
	}

	inline double round_up(double const& d, unsigned int const decimal_place){
		return ceil(d*decimal_place)/decimal_place;
	}
	/*}*/

	/*sign*/
	/*{*/
	template<typename Type> inline constexpr
		int sign_unsigned(Type x){ return Type(0)<x; }

	template<typename Type> inline constexpr
		int sign_signed(Type x){ return (Type(0)<x)-(x<Type(0)); }

	template<typename Type> inline constexpr
		int sign(Type x) { return std::is_signed<Type>()?sign_signed(x):sign_unsigned(x); }
	/*}*/

	/*double real(T)*/
	/*{*/
	inline double real(double const& x){ return x; }

	inline double real(std::complex<double> const& x){ return std::real(x); }
	/*}*/

	/*double imag(T)*/
	/*{*/
	inline double imag(double const& x){ (void)(x); return 0; }

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

	inline void display_progress(double const& step, double const& total, std::string const& msg=" "){
		std::cout<<msg<<100.*step/total<<"%       \r"<<std::flush;
	}

	/*!checks wether a point in inside a polygon (see https://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html) */
	inline bool in_polygon(unsigned int const& nvert, double const* const vertx, double const* const verty, double const& testx, double const& testy){
		unsigned int i(0);
		unsigned int j(0);
		bool c(false);
		for (i = 0, j = nvert-1; i < nvert; j = i++){
			if ( ((verty[i]>testy) != (verty[j]>testy)) && (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) ){ c = !c; }
		}
		return c;
	}

	/*determines if two segments intersect (see from http://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect)*/
	/*{*/
	/*!given three colinear points p, q, r, checks if point q lies between 'pr'*/
	inline bool on_segment(double const* p, double const* q, double const* r){
		return (q[0] <= std::max(p[0], r[0]) && q[0] >= std::min(p[0], r[0]) && q[1] <= std::max(p[1], r[1]) && q[1] >= std::min(p[1], r[1]));
	}

	/*{*//*!Find orientation of ordered triplet (p, q, r).
		   The function returns following values
		   0 --> p, q and r are colinear
		   1 --> Clockwise
		   2 --> Counterclockwise
		   See http://www.geeksforgeeks.org/orientation-3-ordered-points/ // for
		   details of below formula.  *//*}*/
	inline int orientation(double const* p, double const* q, double const* r){
		double val((q[1]-p[1])*(r[0]-q[0]) - (q[0]-p[0])*(r[1]-q[1]));
		if(my::are_equal(val,0)){ return 0; }
		return(val>0)?1:2;
	}

	/*!returns true if line segment 'p1q1' and 'p2q2' intersect*/
	inline bool intersect(double const* p1, double const* q1, double const* p2, double const* q2){
		int o1(my::orientation(p1,q1,p2));
		int o2(my::orientation(p1,q1,q2));
		int o3(my::orientation(p2,q2,p1));
		int o4(my::orientation(p2,q2,q1));
		//General case
		if(o1 != o2 && o3 != o4){ return true; }
		//p1, q1 and p2 are colinear and p2 lies on segment p1q1
		if(o1 == 0 && my::on_segment(p1,p2,q1)){ return true; }
		if(o2 == 0 && my::on_segment(p1,q2,q1)){ return true; }
		//p2, q2 and p1 are colinear and p1 lies on segment p2q2
		if(o3 == 0 && my::on_segment(p2,p1,q2)){ return true; }
		if(o4 == 0 && my::on_segment(p2,q1,q2)){ return true; }
		return false;
	}
	/*}*/

	/*!transforms a double into a fration*/
	inline unsigned int to_fraction(double const& x0, unsigned long long& num, unsigned long long& den, double& sign, double const& err = 1e-15){
		sign = my::sign(x0);
		double g(std::abs(x0));
		unsigned long long a(0);
		unsigned long long b(1);
		unsigned long long c(1);
		unsigned long long d(0);
		unsigned long long s;
		unsigned int iter(0);
		do {
			s = std::floor(g);
			num = a + s*c;
			den = b + s*d;
			a = c;
			b = d;
			c = num;
			d = den;
			g = 1.0/(g-s);
			if(err>std::abs(sign*num/den-x0)){ return iter; }
		} while(iter++<1e6);
		std::cerr<<__PRETTY_FUNCTION__<<" : failed to find a fraction for "<<x0<<std::endl;
		return 0;
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

#ifdef MY_BLAS_ZDOTU
#include <immintrin.h>
#include <emmintrin.h>
	/*!need to redefine this method because for some reason, with an intel
	 * compiler (icpc) the blas::zdotu_ doesn't work properly (at least with
	 * mkl). with gcc, it seems to work fine.*/
	extern "C" {
		inline std::complex<double> zdotu_(unsigned int const& N, std::complex<double> const* const dx, unsigned int ix, std::complex<double> const* const dy, unsigned int iy){
			double const* x = reinterpret_cast<double const*>(dx);
			double const* y = reinterpret_cast<double const*>(dy);

			__m256d dmm0;
			__m256d dmm1;
			__m256d dmm4(_mm256_setzero_pd());
			__m256d dmm5(_mm256_setzero_pd());
			if(ix==1 && iy==1){
				for(unsigned int i(0);i<2*N-2;i+=4){
					dmm0 = _mm256_loadu_pd(x + i);
					dmm1 = _mm256_loadu_pd(y + i);

					/*{would work on Haswell processors
					  dmm4 = _mm256_fmadd_pd(dmm1, dmm0, dmm4);
					  dmm2 = _mm256_permute_pd(dmm1, 0x5);
					  dmm5 = _mm256_fmadd_pd(dmm2, dmm0, dmm5);
					  }*/

					dmm4 = _mm256_add_pd(_mm256_mul_pd(dmm1,dmm0), dmm4);
					dmm1 = _mm256_permute_pd(dmm1, 0x5);
					dmm5 = _mm256_add_pd(_mm256_mul_pd(dmm1,dmm0), dmm5);
				}
			} else {
				for(unsigned int i(0);i<2*N-2;i+=4){
					dmm0 = _mm256_set_pd(x[ix*(i+2)+1], x[ix*(i+2)], x[ix*i+1], x[ix*i]);
					dmm1 = _mm256_set_pd(y[iy*(i+2)+1], y[iy*(i+2)], y[iy*i+1], y[iy*i]);

					/*{would work on Haswell processors
					  dmm4 = _mm256_fmadd_pd(dmm1, dmm0, dmm4);
					  dmm2 = _mm256_permute_pd(dmm1, 0x5);
					  dmm5 = _mm256_fmadd_pd(dmm2, dmm0, dmm5);
					  }*/

					dmm4 = _mm256_add_pd(_mm256_mul_pd(dmm1,dmm0), dmm4);
					dmm1 = _mm256_permute_pd(dmm1, 0x5);
					dmm5 = _mm256_add_pd(_mm256_mul_pd(dmm1,dmm0), dmm5);
				}
			}
			double* re = (double*)&dmm4;
			double* im = (double*)&dmm5;
			std::complex<double> out(re[0]-re[1]+re[2]-re[3],im[0]+im[1]+im[2]+im[3]);
			if(N%2){ out += dx[(N-1)*ix]*dy[(N-1)*iy]; }
			return out;
		}
	}
#else
	extern "C" std::complex<double> zdotu_(unsigned int const& N, std::complex<double> const* const dx, unsigned int const& ix, std::complex<double> const* const dy, unsigned int const& iy);
#endif
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
