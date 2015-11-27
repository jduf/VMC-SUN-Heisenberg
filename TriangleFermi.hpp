#ifndef DEF_TRIANGLEFERMI
#define DEF_TRIANGLEFERMI

#include "Triangle.hpp"

class TriangleFermi: public Triangle<double>{
	public:
		TriangleFermi(System const& s);
		~TriangleFermi() = default;

		void create();
		void check();

	protected:
		void compute_H();
		void display_results();

		Matrix<double> set_ab();
		unsigned int match_pos_in_ab(Vector<double> const& x) const { (void)(x); return 0; }
};
#endif
