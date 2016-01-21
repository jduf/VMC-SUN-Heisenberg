#ifndef DEF_TRIANGLEFREE
#define DEF_TRIANGLEFREE

#include "Triangle.hpp"

class TriangleFree: public Triangle<double>{
	public:
		TriangleFree(System const& s, Vector<double> const& t, Vector<double> const& mu);
		~TriangleFree() = default;

		void create();
		void check();

	protected:
		Vector<double> const t_; //!< hopping term
		Vector<double> const mu_;//!< chemical potential

		void compute_H();
		void display_results();

		Matrix<double> set_ab() const;
		unsigned int match_pos_in_ab(Vector<double> const& x) const;
};
#endif
