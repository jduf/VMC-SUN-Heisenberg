#ifndef DEF_HONEYCOMBFREE
#define DEF_HONEYCOMBFREE

#include "Honeycomb.hpp"

class HoneycombFree: public Honeycomb<double>{
	public:
		HoneycombFree(System const& s, Vector<double> const& t);
		~HoneycombFree() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		Vector<double> const t_; //!< hopping terms

		void compute_H();
		void display_results();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;
};
#endif
