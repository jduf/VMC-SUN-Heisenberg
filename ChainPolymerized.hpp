#ifndef DEF_CHAINPOLYMERIZED
#define DEF_CHAINPOLYMERIZED

#include "Chain.hpp"

class ChainPolymerized: public Chain<double> {
	public:
		ChainPolymerized(unsigned int N, unsigned int n, unsigned int m);
		~ChainPolymerized();

		void create(double delta);
		void study();
		void save(Write& w);

	private:
		void compute_P(Matrix<double>& P);
		void compute_T();

		double delta_;
};
#endif
