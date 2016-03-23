#ifndef DEF_SQUARET4X2
#define DEF_SQUARET4X2

#include "Square.hpp"

class SquareT4x2: public Square<double>{
	public:
		SquareT4x2(System const& s, Vector<double> const& t);
		~SquareT4x2() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	protected:
		Vector<double> const t_; //!< hopping terms

		void compute_H();
		void display_results();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;
};
#endif
