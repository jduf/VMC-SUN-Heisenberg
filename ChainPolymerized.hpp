#ifndef DEF_CHAINPOLYMERIZED
#define DEF_CHAINPOLYMERIZED

#include "Chain.hpp"

/*!{Creates a chain with a weaker hopping parameter every N/m link*/
/*
 * To properly solve the degeneracy problem, this wavefunction selects the
 * eigenvector |0>. This works for all cases but when N/m and nm/N are even.
 */
/*}*/
class ChainPolymerized: public Chain<double> {
	public:
		ChainPolymerized(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc, double delta);
		~ChainPolymerized(){}

		void create();
		void save() const;
		void check();
		
	private:
		double delta_;//!< polymerization parameter 

		void compute_H();
		std::string extract_level_7();
		std::string extract_level_6();
};
#endif
