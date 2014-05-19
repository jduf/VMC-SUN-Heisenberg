#ifndef DEF_SQUAREPIFLUX
#define DEF_SQUAREPIFLUX

#include "Square.hpp"

class SquarePiFlux: public Square<std::complex<double> >{
	public:
		SquarePiFlux(unsigned int N, unsigned int n, unsigned int m, int bc);
		~SquarePiFlux();

		bool create(double x);
		void check();

	protected:
		void compute_T();
};
#endif
