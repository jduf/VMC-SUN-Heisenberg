#ifndef DEF_SQUAREFREEFLUX
#define DEF_SQUAREFREEFLUX

#include "Square.hpp"

class SquareFreeFlux: public Square<std::complex<double> >{
	public:
		SquareFreeFlux(System const& s, Vector<double> const& phi);
		~SquareFreeFlux() = default;

		void create();
		void check();

	protected:
		Vector<double> const phi_;

		void compute_H();
		void display_results();

		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
