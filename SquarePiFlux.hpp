#ifndef DEF_SQUAREPIFLUX
#define DEF_SQUAREPIFLUX

#include "Square.hpp"

class SquarePiFlux: public Square<std::complex<double> >{
	public:
		SquarePiFlux(Parseur& P);
		~SquarePiFlux();

	protected:
		void compute_T();
		void save();
};
#endif

