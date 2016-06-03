#ifndef DEF_HONEYCOMBP00
#define DEF_HONEYCOMBP00

#include "Honeycomb.hpp"

/*{*//*!Plaquette wavefunction with pi00-flux configuration if td<0 (each
	   pi-flux hexagon is surrounded by 0-flux hexagons), pipipi-flux
	   otherwise.
	   *//*}*/
class Honeycombp00: public Honeycomb<double>{
	public:
		Honeycombp00(System const& s, double const& td);
		~Honeycombp00() = default;

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

		std::string get_mc_run_command() const;
};
#endif
