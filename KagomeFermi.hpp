#ifndef DEF_KAGOMEFERMI
#define DEF_KAGOMEFERMI

#include "Kagome.hpp"

class KagomeFermi: public Kagome<double>{
	public:
		KagomeFermi(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc);
		~KagomeFermi(){}

		void create();
		void check();

	protected:
		void compute_T();
		void lattice();
};
#endif
