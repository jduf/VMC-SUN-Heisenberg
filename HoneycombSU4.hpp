#ifndef DEF_HONEYCOMBSU4
#define DEF_HONEYCOMBSU4

#include "Honeycomb.hpp"

class HoneycombSU4: public Honeycomb<double>{
	public:
		HoneycombSU4(System const& s);
		~HoneycombSU4() = default;

		void create();
		void check();

	protected:
		void compute_H();
		void display_results();

		Matrix<double> set_ab();
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif

