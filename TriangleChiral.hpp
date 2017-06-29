#ifndef DEF_TRIANGLECHIRAL
#define DEF_TRIANGLECHIRAL

#include "Triangle.hpp"

class TriangleChiral: public Triangle<std::complex<double> >{
	public:
		TriangleChiral(System const& s, double const& phi);
		~TriangleChiral() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		double const phi_; //!< flux per triangular plaquette

		void compute_H();
		void display_results();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab(unsigned int const& k) const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;

		std::string extract_level_2();
};
#endif
