#ifndef DEF_TRIANGLEPHI
#define DEF_TRIANGLEPHI

#include "Triangle.hpp"

class TrianglePhi: public Triangle<std::complex<double> >{
	public:
		TrianglePhi(unsigned int N, unsigned int n, unsigned int m);
		~TrianglePhi();

		void create(double phi);
		void check();
		void study();
		void save(Write& w) const;

	protected:
		double phi_;
		
		void compute_T();
		void compute_P(Matrix<std::complex<double> >& Px, Matrix<std::complex<double> >& Py);
};
#endif

