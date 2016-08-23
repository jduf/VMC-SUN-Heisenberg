#ifndef DEF_TRIANGLET3X2
#define DEF_TRIANGLET3X2

#include "Triangle.hpp"

class TriangleT3x2: public Triangle<double>{
	public:
		TriangleT3x2(System const& s, Vector<double> const& t);
		~TriangleT3x2() = default;

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
