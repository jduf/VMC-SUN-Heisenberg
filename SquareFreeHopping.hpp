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
		void display_results();
		Vector<double> const t_;

		/*!Set the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
