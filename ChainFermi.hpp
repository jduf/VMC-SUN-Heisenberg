#ifndef DEF_CHAINFERMI
#define DEF_CHAINFERMI

#include "Chain.hpp"

class ChainFermi: public Chain<double>{
	public:
		ChainFermi(Container const& P);
		~ChainFermi();

		void create(double x=0);
		void study();
		void save();

	private:
		void compute_P(Matrix<double>& P);
		void compute_T();
};
#endif
