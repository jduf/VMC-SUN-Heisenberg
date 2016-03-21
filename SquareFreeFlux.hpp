#ifndef DEF_SQUAREFREEFLUX
#define DEF_SQUAREFREEFLUX

#include "Square.hpp"

class SquareFreeFlux: public Square<std::complex<double> >{
	public:
		SquareFreeFlux(System const& s, Vector<double> const& t, Vector<double> const& phi);
		~SquareFreeFlux() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	protected:
		Vector<double> const t_; //!< hopping terms
		Vector<double> const phi_;

		void compute_H();
		void display_results();

		/*!Set the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;
};
#endif
