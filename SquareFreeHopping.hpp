#ifndef DEF_SQUAREFREEREAL
#define DEF_SQUAREFREEREAL

#include "Square.hpp"

class SquareFreeHopping: public Square<double>{
	public:
		SquareFreeHopping(System const& s, Vector<double> const& t);
		~SquareFreeHopping() = default;

		void create();
		void check();

	protected:
		void compute_H();
		void lattice(std::string const& path, std::string const& filename);
		Vector<double> const t_;

		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
