#ifndef DEF_TRIANGLEFERMI
#define DEF_TRIANGLEFERMI

#include "Triangle.hpp"

class TriangleFermi: public Triangle<double>{
	public:
		TriangleFermi(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc);
		~TriangleFermi(){}

		void create();
		void check();

	protected:
		void compute_H();
		void lattice();

		unsigned int match_pos_in_ab(Vector<double> const& x) const;
		Matrix<double> set_ab();
};
#endif
