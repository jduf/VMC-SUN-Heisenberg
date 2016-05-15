#ifndef DEF_SQUARECHIRAL
#define DEF_SQUARECHIRAL

#include "Square.hpp"

class SquareChiral: public Square<std::complex<double> >{
	public:
		SquareChiral(System const& s, double const& phi);
		~SquareChiral() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		double const phi_; //!< flux per square plaquette

		void compute_H();
		void display_results();

		/*!Set the unit cell's vectors*/
		Matrix<double> set_ab(unsigned int const& ref3, unsigned int const& k) const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;
};
#endif
