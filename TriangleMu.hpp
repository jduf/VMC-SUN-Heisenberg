#ifndef DEF_TRIANGLEMU
#define DEF_TRIANGLEMU

#include "Triangle.hpp"

class TriangleMu: public Triangle<double>{
	public:
		TriangleMu(System const& s, double const& mu);
		~TriangleMu() = default;

		void create();
		void check();

	protected:
		double const mu_;//!< chemical potential

		void compute_H(unsigned int const& c);
		void display_results();

		Matrix<double> set_ab() const;
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
