#ifndef DEF_CHAINFERMI
#define DEF_CHAINFERMI

#include "Chain.hpp"

class ChainFermi: public Chain<double>{
	public:
		ChainFermi(unsigned int N, unsigned int n, unsigned int m, int bc);
		~ChainFermi();

		bool create(double x=0);
		void check();
		void study(double E, double DeltaE, Vector<double> const& corr, std::string save_in);

	private:
		void compute_P(Matrix<double>& P);
		void compute_T();
};
#endif
