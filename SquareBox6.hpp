#ifndef DEF_SQUAREBOX6
#define DEF_SQUAREBOX6

#include "Square.hpp"

class SquareBox6: public Square<std::complex<double> >{
	public:
		SquareBox6(System const& s, Vector<double> const& t);
		~SquareBox6() = default;

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
