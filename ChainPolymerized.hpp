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
		ChainPolymerized(System const& s, Vector<double> const& t);
		~ChainPolymerized() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		Vector<double> t_;//!< polymerization parameter

		void compute_H();
		std::string extract_level_8();
		std::string extract_level_7();

		void energy_bound();
		unsigned int set_spuc(Vector<double> const& t, unsigned int const& spuc);
};
#endif
