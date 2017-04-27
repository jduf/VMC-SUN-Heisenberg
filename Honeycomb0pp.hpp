#ifndef DEF_HONEYCOMB0PP
#define DEF_HONEYCOMB0PP

#include "Honeycomb.hpp"

/*{*//*!Plaquette wavefunction with 0pipi-flux configuration if td<0 (each
	   0-flux hexagon is surrounded by pi-flux hexagons), 000-flux otherwise.
	   *//*}*/
class Honeycomb0pp: public Honeycomb<double>{
	public:
		Honeycomb0pp(System const& s, double const& td);
		~Honeycomb0pp() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		double const td_;

		void compute_H();
		void display_results();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;
};
#endif
