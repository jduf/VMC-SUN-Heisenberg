#ifndef DEF_SQUAREPIFLUX
#define DEF_SQUAREPIFLUX

#include "Square.hpp"

#include <complex>

class SquarePiFlux: public Square<std::complex<double> >{
	public:
		SquarePiFlux(Parseur& P);
		~SquarePiFlux();

	protected:
		void compute_T();
		void save(std::string filename);
};

#endif

