#ifndef DEF_SQUAREACSL
#define DEF_SQUAREACSL

#include "Square.hpp"

class SquareACSL: public Square<std::complex<double> >{
	public:
		SquareACSL(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc, Vector<double> const& t);
		~SquareACSL() = default;

		void create();
		void check();

	protected:
		void compute_H();
		void lattice();
		Vector<double> const t_;

		Matrix<double> set_ab(unsigned int const& spuc);
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
