#ifndef DEF_TRIANGLECHIRAL
#define DEF_TRIANGLECHIRAL

#include "Triangle.hpp"

class TriangleChiral: public Triangle<std::complex<double> >{
	public:
		TriangleChiral(System const& s);
		~TriangleChiral() = default;

		void create();
		void check();

	protected:
		double const phi_;

		void compute_H();
		void display_results();

		Matrix<double> set_ab() const;
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
