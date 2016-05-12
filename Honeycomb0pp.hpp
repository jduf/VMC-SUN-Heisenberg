#ifndef DEF_HONEYCOMB0PP
#define DEF_HONEYCOMB0PP

#include "Honeycomb.hpp"

/*{*//*!Plaquette wavefunction with possible 0-pi-pi flux configuration
	   Three different flux configurations (fc) are possible : taking links
	   0->1, 2->3, 4->5 as reference and closing each hexagon with the same
	   orientation, one has :

	   | fc   : 0  |  1 | 2
	   | 0->1 : pi | pi | 0
	   | 2->3 : pi | 0  | pi
	   | 4->5 : 0  | pi | pi

	   They are all equivalent but the overlap matrix after projection is not
	   equal to the identity.
	   *//*}*/
class Honeycomb0pp: public Honeycomb<double>{
	public:
		Honeycomb0pp(System const& s, double const& td, unsigned int const& fc);
		~Honeycomb0pp() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		double const td_;
		unsigned int const fc_;

		void compute_H();

		void display_results();
		void lattice();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;

		std::string extract_level_2();
};
#endif
