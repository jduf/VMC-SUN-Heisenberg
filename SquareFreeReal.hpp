#ifndef DEF_SQUAREFREEREAL
#define DEF_SQUAREFREEREAL

#include "Square.hpp"

class SquareFreeReal: public Square<double>{
	public:
		SquareFreeReal(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc, Vector<double> const& t, Vector<double> const& mu);
		~SquareFreeReal() = default;

		void create();
		void check();

	protected:
		void compute_H(unsigned int const& c);
		void lattice();
		Vector<double> const t_;
		Vector<double> const mu_;

		unsigned int match_pos_in_ab(Vector<double> const& x) const;
		Matrix<double> set_ab();
};
#endif
