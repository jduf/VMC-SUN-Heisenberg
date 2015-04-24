#ifndef DEF_KAGOMEVBC
#define DEF_KAGOMEVBC

#include "Kagome.hpp"

class KagomeVBC: public Kagome<std::complex<double> >{
	public:
		KagomeVBC(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc);
		~KagomeVBC() = default;

		void create();
		void check();

	protected:
		void compute_H();
		void lattice();
		std::string extract_level_7();
		std::string extract_level_6();
		std::string extract_level_4();
};
#endif
