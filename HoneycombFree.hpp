#ifndef DEF_HONEYCOMBFREE
#define DEF_HONEYCOMBFREE

#include "Honeycomb.hpp"

class HoneycombFree: public Honeycomb<double>{
	public:
		HoneycombFree(System const& s, Vector<double> t);
		~HoneycombFree() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	protected:
		Vector<double> t_;

		void compute_H();

		void display_results();
		void lattice();

		Matrix<double> set_ab() const;
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
