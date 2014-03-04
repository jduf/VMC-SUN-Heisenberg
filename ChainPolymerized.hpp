#ifndef DEF_CHAINPOLYMERIZED
#define DEF_CHAINPOLYMERIZED

#include "Chain.hpp"

class ChainPolymerized: public Chain<double> {
	public:
		ChainPolymerized(unsigned int N, unsigned int n, unsigned int m, int bc);
		~ChainPolymerized();

		unsigned int create(double delta);
		void check();
		void study();
		void save(Write& w) const;

	private:
		void compute_P(Matrix<double>& P);
		void compute_T();

		double delta_;
};
#endif
