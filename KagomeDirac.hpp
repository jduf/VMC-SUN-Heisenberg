#ifndef DEF_KAGOMEDIRAC
#define DEF_KAGOMEDIRAC

#include "Kagome.hpp"

class KagomeDirac: public Kagome<double>{
	public:
		KagomeDirac(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc);
		~KagomeDirac();

		void create();
		void check();

	protected:
		void compute_T();
		void compute_P(Matrix<double>& Px, Matrix<double>& Py);
		void lattice();
};
#endif
