#ifndef DEF_HONEYCOMBPIFLUX
#define DEF_HONEYCOMBPIFLUX

#include "Honeycomb.hpp"

class HoneycombPiFlux: public Honeycomb<double>{
	public:
		HoneycombPiFlux(System const& s);
		~HoneycombPiFlux() = default;

		void create();
		void check();

	protected:
		void compute_H();
		void display_results();

		Matrix<double> set_ab();
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif

