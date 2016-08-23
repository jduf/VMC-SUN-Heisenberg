#ifndef DEF_SQUARET4X4
#define DEF_SQUARET4X4

#include "Square.hpp"

class SquareT4x4: public Square<double>{
	public:
		SquareT4x4(System const& s, Vector<double> const& t);
		~SquareT4x4() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		Vector<double> const t_; //!< hopping terms

		void compute_H();
		void display_results();

		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;
};
#endif
