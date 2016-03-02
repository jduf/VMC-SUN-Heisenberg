#ifndef DEF_SQUAREFERMI
#define DEF_SQUAREFERMI

#include "Square.hpp"

class SquareFermi: public Square<double>{
	public:
		SquareFermi(System const& s);
		~SquareFermi() = default;

		void create();
		void check();

	protected:
		void compute_H();
		void display_results();
		void lattice();

		/*!Set the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const { (void)(x); return 0; }
};
#endif
