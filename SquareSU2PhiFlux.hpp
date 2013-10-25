#ifndef DEF_SQUARESU2PHIFLUX
#define DEF_SQUARESU2PHIFLUX

#include "Square.hpp"

class SquareSU2PhiFlux: public Square<std::complex<double> >{
	public:
		SquareSU2PhiFlux(Parseur& P);
		~SquareSU2PhiFlux();

	protected:
		double const phi_;

		void compute_T();
		void save();

		void compute_P();
		void band_structure();
		void lattice();
};
#endif

