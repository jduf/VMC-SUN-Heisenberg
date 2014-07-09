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
		std::string extract_level_7();
		std::string extract_level_6();
		std::string extract_level_4();
		std::string extract_level_3();
};
#endif
