#ifndef DEF_CHAINPOLYMERIZED
#define DEF_CHAINPOLYMERIZED

#include "Chain.hpp"

class ChainPolymerized: public Chain<double> {
	public:
		ChainPolymerized(unsigned int N, unsigned int n, unsigned int m, int bc);
		~ChainPolymerized();

		bool create(double delta);
		void check();
		void save(Write& w) const;
		void study(double E, double DeltaE, Vector<double> const& corr, std::string save_in);

	private:
		void compute_P(Matrix<double>& P);
		void compute_T();

		double delta_;
};
#endif
