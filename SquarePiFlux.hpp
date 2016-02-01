#ifndef DEF_SQUAREPIFLUX
#define DEF_SQUAREPIFLUX

#include "Square.hpp"

class SquarePiFlux: public Square<std::complex<double> >{
	public:
		SquarePiFlux(System const& s);
		~SquarePiFlux() = default;

		void create();
		void check();

	protected:
		void compute_H();
		void display_results();

		std::string extract_level_7();
		std::string extract_level_3();

		/*!Set the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
