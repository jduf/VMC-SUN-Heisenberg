#ifndef DEF_TRIANGLEPLAQUETTE
#define DEF_TRIANGLEPLAQUETTE

#include "Triangle.hpp"

class TrianglePlaquette: public Triangle<double>{
	public:
		TrianglePlaquette(System const& s, double const& t);
		~TrianglePlaquette() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	protected:
		double const t_; //!< hopping term

		void compute_H();
		void display_results();

		Matrix<double> set_ab() const;
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
