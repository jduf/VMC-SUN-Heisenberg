#ifndef DEF_SQUARESU2PHIFLUX
#define DEF_SQUARESU2PHIFLUX

#include "Square.hpp"

class SquareSU2PhiFlux: public Square<std::complex<double> >{
	public:
		SquareSU2PhiFlux(Container const& param);
		~SquareSU2PhiFlux();

		void create(double phi);
		void study();
		void save();
		void get_param(Container& param);

	protected:
		double phi_;

		void compute_T();
		void compute_P(Matrix<std::complex<double> >& Px, Matrix<std::complex<double> >& Py);
		void lattice();
};
#endif

