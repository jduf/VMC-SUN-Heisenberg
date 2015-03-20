#ifndef DEF_SQUAREJASTROW
#define DEF_SQUAREJASTROW

#include "Square.hpp"

class SquareJastrow: public Square<double>{
	public:
		SquareJastrow(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc, Matrix<double> const& nu);
		~SquareJastrow(){}

		void create();
		void check();
		void save() const;

	protected:
		void compute_nn();
		void compute_sublattice();
		void compute_omega_cc();
		void lattice();

		unsigned int match_pos_in_ab(Vector<double> const& x) const;
		Matrix<double> set_ab();
};
#endif

