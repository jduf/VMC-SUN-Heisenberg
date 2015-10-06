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
		void lattice(std::string const& path, std::string const& filename);

		std::string extract_level_7();
		std::string extract_level_3();

		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
