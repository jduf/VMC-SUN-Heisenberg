#ifndef DEF_HONEYCOMBCHIRAL
#define DEF_HONEYCOMBCHIRAL

#include "Honeycomb.hpp"

class HoneycombChiral: public Honeycomb<std::complex<double> >{
	public:
		HoneycombChiral(System const& s, double const& phi);
		~HoneycombChiral() = default;

		void create();
		void check();

	protected:
		double const phi_; //!< flux per hexagonal plaquette

		void compute_H();
		void display_results();
		
		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
