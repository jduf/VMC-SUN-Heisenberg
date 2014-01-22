#ifndef DEF_TRIANGLEPHI
#define DEF_TRIANGLEPHI

#include "Triangle.hpp"

class TrianglePhi: public Triangle<std::complex<double> >{
	public:
		TrianglePhi(Parseur& P);
		~TrianglePhi();

		void save();
	protected:
		double phi_;
		
		void compute_P(Matrix<std::complex<double> >& Px, Matrix<std::complex<double> >& Py);
		void compute_T();
};
#endif

