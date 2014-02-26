#ifndef DEF_CHAINDIMERIZED
#define DEF_CHAINDIMERIZED

#include "Chain.hpp"

class ChainDimerized: public Chain<double> {
	public:
		ChainDimerized(unsigned int N, unsigned int n, unsigned int m);
		~ChainDimerized();

		void create(double delta);
		void study();
		void save();

	private:
		void compute_P(Matrix<double>& P);
		void compute_T();

		double delta_;
};
#endif
