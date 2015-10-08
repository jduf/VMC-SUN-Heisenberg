#ifndef DEF_SQUAREFREEFLUX
#define DEF_SQUAREFREEFLUX

#include "Square.hpp"

class SquareFreeFlux: public Square<std::complex<double> >{
	public:
		SquareFreeFlux(System const& s, Vector<double> const& phi);
		~SquareFreeFlux() = default;

		void create(unsigned int const& which_observables);
		void check();

	protected:
		void compute_H();
		void lattice(std::string const& path, std::string const& filename);
		Vector<double> const phi_;

		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
