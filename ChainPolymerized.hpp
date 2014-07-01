#ifndef DEF_CHAINPOLYMERIZED
#define DEF_CHAINPOLYMERIZED

#include "Chain.hpp"

class ChainPolymerized: public Chain<double> {
	public:
		ChainPolymerized(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc, double delta);
		~ChainPolymerized();

		void create();
		void check();
		void save(IOFiles& w) const;
		
	private:
		void compute_P(Matrix<double>& P);
		void compute_T();

		double delta_;
};
#endif
