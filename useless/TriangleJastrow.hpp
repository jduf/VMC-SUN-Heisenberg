#ifndef DEF_TRIANGLEJASTROW
#define DEF_TRIANGLEJASTROW

#include "Triangle.hpp"

class TriangleJastrow: public Triangle<double>{
	public:
		TriangleJastrow(System const& s, Matrix<double> const& nu);
		~TriangleJastrow() = default;

		void create();
		void check();
		void save_param(IOFiles& w) const;

	protected:
		void display_results();

		/*!Set the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const;
};
#endif

