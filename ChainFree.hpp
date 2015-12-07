#ifndef DEF_CHAINFREE
#define DEF_CHAINFREE

#include "Chain.hpp"

/*!{Creates a chain with a weaker hopping parameter every N/m link*/
/*
 * To properly solve the degeneracy problem, this wavefunction selects the
 * eigenvector |0>. This works for all cases but when N/m and nm/N are even.
 */
/*}*/
class ChainFree: public Chain<double> {
	public:
		ChainFree(System const& s, Vector<double> const& t, Vector<double> const& mu);
		~ChainFree() = default;

		void create();
		void save_param(IOFiles& w) const;
		void check();

	private:
		Vector<double> t_; //!< polymerization parameter
		Vector<double> mu_;//!< chemical potential

		void compute_H();
		unsigned int set_spuc(Vector<double> const& t, Vector<double> const& mu, unsigned int const& spuc);

		void display_results();
		void energy_bound();
};
#endif
