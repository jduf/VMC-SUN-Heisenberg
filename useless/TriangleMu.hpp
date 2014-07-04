#ifndef DEF_TRIANGLEMU
#define DEF_TRIANGLEMU

#include "Triangle.hpp"

class TriangleMu: public Triangle<double>{
	public:
		TriangleMu(unsigned int N, unsigned int n, unsigned int m);
		~TriangleMu();

		void create(double mu);
		void check();
		void study();
		void save(Write& w) const;

	protected:
		double mu_;

		void compute_T(unsigned int alpha);
		void compute_P(Matrix<double>& Px, Matrix<double>& Py);
		void lattice();
};
#endif

