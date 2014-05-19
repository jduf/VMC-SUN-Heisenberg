#ifndef DEF_SQUAREFERMI
#define DEF_SQUAREFERMI

#include "Square.hpp"

class SquareFermi: public Square<double>{
	public:
		SquareFermi(unsigned int N, unsigned int n, unsigned int m);
		~SquareFermi();

		void create(double x);
		void check();

	protected:
		void compute_T();
		void compute_P(Matrix<double>& Px, Matrix<double>& Py);
		void lattice();
};
#endif
