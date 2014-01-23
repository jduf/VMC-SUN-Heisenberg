#ifndef DEF_CHAINFERMI
#define DEF_CHAINFERMI

#include "Chain.hpp"

class ChainFermi: public Chain<double>{
	public:
		ChainFermi(Parseur& P);
		~ChainFermi();

		void save();

	private:
		void compute_P(Matrix<double>& P);
		void compute_T();
};
#endif
