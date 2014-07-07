#ifndef DEF_CHAINPOLYMERIZED
#define DEF_CHAINPOLYMERIZED

#include "Chain.hpp"

class ChainPolymerized: public Chain<double> {
	public:
		ChainPolymerized(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc, double delta);
		~ChainPolymerized(){}

		void create();
		void save(IOFiles& w) const;
		void check();
		std::string analyse(IOSystem const& t);
		
	private:
		double delta_;

		void compute_T();
};
#endif
