#ifndef DEF_SQUARECHIRAL
#define DEF_SQUARECHIRAL

#include "Square.hpp"

class SquareChiral: public Square<std::complex<double> >{
	public:
		SquareChiral(System const& s);
		~SquareChiral() = default;

		void create();
		void check();

	protected:
		double const phi_;

		void compute_H();
		void display_results();

		/*!Set the unit cell's vectors*/
		Matrix<double> set_ab(unsigned int const& k) const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
