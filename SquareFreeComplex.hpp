#ifndef DEF_SQUAREFREECOMPLEX
#define DEF_SQUAREFREECOMPLEX

#include "Square.hpp"

class SquareFreeComplex: public Square<std::complex<double> >{
	public:
		SquareFreeComplex(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc, Vector<double> const& t, Vector<double> const& mu, Vector<double> const& phi);
		~SquareFreeComplex(){}

		void create();
		void check();

	protected:
		void compute_H(unsigned int const& c);
		void lattice();
		Vector<double> const t_;
		Vector<double> const mu_;
		Vector<double> const phi_;

		Matrix<double> set_ab();
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
