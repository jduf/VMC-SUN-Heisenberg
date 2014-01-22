#ifndef DEF_SQUARESU2PHIFLUX
#define DEF_SQUARESU2PHIFLUX

#include "Square.hpp"

class SquareSU2PhiFlux: public Square<std::complex<double> >{
	public:
		SquareSU2PhiFlux(Parseur& P);
		~SquareSU2PhiFlux();

		void save();

	protected:
		double const phi_;

		void compute_T();
		void compute_P(Matrix<std::complex<double> >& Px, Matrix<std::complex<double> >& Py);
		void lattice();
};
#endif

