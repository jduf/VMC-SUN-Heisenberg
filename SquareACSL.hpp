#ifndef DEF_SQUAREACSL
#define DEF_SQUAREACSL

#include "Square.hpp"

class SquareACSL: public Square<std::complex<double> >{
	public:
		SquareACSL(System const& s, Vector<double> const& t);
		~SquareACSL() = default;

		void create();
		void check();

	protected:
		void compute_H();
		void lattice(std::string const& path);
		Vector<double> const t_;

		Matrix<double> set_ab(unsigned int const& spuc);
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
