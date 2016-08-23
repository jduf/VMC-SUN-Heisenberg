#ifndef DEF_TRIANGLEPHI
#define DEF_TRIANGLEPHI

#include "Triangle.hpp"

class TrianglePhi: public Triangle<std::complex<double> >{
	public:
		TrianglePhi(System const& s, double const& phi);
		~TrianglePhi() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		double const phi_;

		void compute_H();
		void display_results();
		
		/*!Sets the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int unit_cell_index(Vector<double> const& x) const { (void)(x); return 0; }
};
#endif
