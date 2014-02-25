#ifndef DEF_SQUAREPIFLUX
#define DEF_SQUAREPIFLUX

#include "Square.hpp"

class SquarePiFlux: public Square<std::complex<double> >{
	public:
		SquarePiFlux(Container const& P);
		~SquarePiFlux();

		void create(double x);
		void study();
		void save();

	protected:
		void compute_T();
};
#endif

