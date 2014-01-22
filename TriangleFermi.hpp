#ifndef DEF_TRIANGLEFERMI
#define DEF_TRIANGLEFERMI

#include "Triangle.hpp"

class TriangleFermi: public Triangle<double>{
	public:
		TriangleFermi(Parseur& P);
		~TriangleFermi();

		void save();

	protected:
		void compute_T();
		void compute_P(Matrix<double>& Px, Matrix<double>& Py);
		void lattice();
};
#endif
