#ifndef DEF_HONEYCOMBPLAQUETTE
#define DEF_HONEYCOMBPLAQUETTE

#include "Honeycomb.hpp"

class HoneycombPlaquette: public Honeycomb<double>{
	public:
		HoneycombPlaquette(System const& s, Vector<double> const& t);
		~HoneycombPlaquette() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	protected:
		Vector<double> const t_; //!< hopping terms

		void compute_H();

		void display_results();
		void lattice();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
