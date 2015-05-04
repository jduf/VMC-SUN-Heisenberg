#ifndef DEF_SQUAREPIFLUX
#define DEF_SQUAREPIFLUX

#include "Square.hpp"

class SquarePiFlux: public Square<std::complex<double> >{
	public:
		SquarePiFlux(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc);
		~SquarePiFlux() = default;

		void create();
		void check();

	protected:
		void compute_H();
		void lattice();

		std::string extract_level_7();
		std::string extract_level_3();

		Matrix<double> set_ab();
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
