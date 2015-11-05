#ifndef DEF_SQUAREJASTROW
#define DEF_SQUAREJASTROW

#include "Square.hpp"

class SquareJastrow: public Square<double>{
	public:
		SquareJastrow(System const& s, Matrix<double> const& nu);
		~SquareJastrow() = default;

		void create();
		void check();
		void save_input(IOFiles& w) const;

	protected:
		void compute_nn();
		void compute_sublattice();
		void compute_omega_cc();
		void display_results();

		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
