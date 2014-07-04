#ifndef DEF_SQUAREMU
#define DEF_SQUAREMU

#include "Square.hpp"

class SquareMu: public Square<double>{
	public:
		SquareMu(unsigned int N, unsigned int n, unsigned int m);
		~SquareMu();

		void create(double mu);
		void check();
		void save(Write& w) const;

	protected:
		double mu_;

		void compute_T(unsigned int alpha);
		void compute_P(Matrix<double>& Px, Matrix<double>& Py);
		void lattice();
};
#endif

