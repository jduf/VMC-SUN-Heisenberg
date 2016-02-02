#ifndef DEF_SQUAREFREEREAL
#define DEF_SQUAREFREEREAL

#include "Square.hpp"

class SquareFreeHopping: public Square<double>{
	public:
		SquareFreeHopping(System const& s, Vector<double> const& t, Vector<double> const& mu);
		~SquareFreeHopping() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	protected:
		Vector<double> const t_; //!< hopping terms
		Vector<double> const mu_;//!< chemical potentials
		std::vector<std::pair<unsigned int,unsigned int> > extra_links_;

		void compute_H(unsigned int const& c);
		void display_results();

		/*!Set the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
