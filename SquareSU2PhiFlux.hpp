#ifndef DEF_SQUARESU2PHIFLUX
#define DEF_SQUARESU2PHIFLUX

#include "Square.hpp"
#include "PSTricks.hpp"

#include <complex>

class SquareSU2PhiFlux: public Square<std::complex<double> >{
	public:
		SquareSU2PhiFlux(Parseur& P);
		~SquareSU2PhiFlux();

	protected:
		double const phi_;
		void compute_T();
		void compute_P();
		void compute_band_structure();
		void show_bound();
		void save(std::string filename);
};
#endif

