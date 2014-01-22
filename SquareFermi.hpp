#ifndef DEF_SQUAREFERMI
#define DEF_SQUAREFERMI

#include "Square.hpp"

class SquareFermi: public Square<double>{
	public:
		SquareFermi(Parseur& P);
		~SquareFermi();

		void save();

	protected:
		void compute_T();
		void compute_P(Matrix<double>& Px, Matrix<double>& Py);
		void lattice();
};
#endif
