#ifndef DEF_SQUAREFREE
#define DEF_SQUAREFREE

#include "Square.hpp"

class SquareFree: public Square<double>{
	public:
		SquareFree(System const& s, Vector<double> const& t, Vector<double> const& mu);
		~SquareFree() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		Vector<double> const t_; //!< hopping terms
		Vector<double> const mu_;//!< chemical potentials

		void init_additional_links();
		void compute_H(unsigned int const& c);
		void display_results();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab(unsigned int const& ref3) const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;
};
#endif
