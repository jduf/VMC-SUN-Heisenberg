#ifndef DEF_CHAINFERMI
#define DEF_CHAINFERMI

#include "Chain.hpp"

class ChainFermi: public Chain<double>{
	public:
		ChainFermi(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc);
		~ChainFermi(){}

		void create();
		void check();

	private:
		void compute_T();
};
#endif
