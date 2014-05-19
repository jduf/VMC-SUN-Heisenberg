#ifndef DEF_SQUARESU2PHIFLUX
#define DEF_SQUARESU2PHIFLUX

#include "Square.hpp"

class SquareSU2PhiFlux: public Square<std::complex<double> >{
	public:
		SquareSU2PhiFlux(unsigned int N, unsigned int n, unsigned int m);
		~SquareSU2PhiFlux();

		void create(double phi);
		void check();
		void save(Write& w) const;

	protected:
		double phi_;

		void compute_T();
		void compute_P(Matrix<std::complex<double> >& Px, Matrix<std::complex<double> >& Py);
		void lattice();
};
#endif

