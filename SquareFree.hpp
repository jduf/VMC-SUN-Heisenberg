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

	protected:
		Vector<double> const t_; //!< hopping terms
		Vector<double> const mu_;//!< chemical potentials

		void init_additional_links();
		void compute_H(unsigned int const& c);
		void display_results();

		/*!Set the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
