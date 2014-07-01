#ifndef DEF_KAGOMEVBC
#define DEF_KAGOMEVBC

#include "Kagome.hpp"

class KagomeVBC: public Kagome<std::complex<double> >{
	public:
		KagomeVBC(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc);
		~KagomeVBC();

		void create();
		void check();

	protected:
		void compute_T();
		void compute_P(Matrix<std::complex<double> >& Px, Matrix<std::complex<double> >& Py);
		void lattice();
};
#endif
