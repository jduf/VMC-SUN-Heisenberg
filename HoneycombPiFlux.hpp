#ifndef DEF_HONEYCOMBPIFLUX
#define DEF_HONEYCOMBPIFLUX

#include "Honeycomb.hpp"

class HoneycombPiFlux: public Honeycomb<double>{
	public:
		HoneycombPiFlux(System const& s);
		~HoneycombPiFlux() = default;

		void create();
		void check();

	private:
		void compute_H();
		void display_results();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;
};
#endif

