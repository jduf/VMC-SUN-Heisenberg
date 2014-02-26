#ifndef DEF_SQUAREPIFLUX
#define DEF_SQUAREPIFLUX

#include "Square.hpp"

class SquarePiFlux: public Square<std::complex<double> >{
	public:
		SquarePiFlux(unsigned int N, unsigned int n, unsigned int m);
		~SquarePiFlux();

		void create(double x);
		void study();

	protected:
		void compute_T();
};
#endif
