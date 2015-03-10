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

		static Matrix<double> set_ab(){
			Matrix<double> tmp(2,2);
			tmp(0,0) = 1;
			tmp(1,0) = 0;
			tmp(0,1) = 0;
			tmp(1,1) = 1;
			return tmp;
		}
};
#endif
