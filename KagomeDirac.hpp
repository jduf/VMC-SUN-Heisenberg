#ifndef DEF_KAGOMEDIRAC
#define DEF_KAGOMEDIRAC

#include "Kagome.hpp"

class KagomeDirac: public Kagome<double>{
	public:
		KagomeDirac(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc);
		~KagomeDirac(){}

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
