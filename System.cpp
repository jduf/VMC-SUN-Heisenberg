#include "System.hpp"

/*constructors and destructor*/
/*{*/
System::System(unsigned int const& N, unsigned int const& n, unsigned int const& m, int const& bc, Vector<unsigned int> const& ref):
	ref_(ref),
	n_(n),
	N_(N), 
	m_(m),
	M_((m_*n_)/N_), 
	bc_(bc)
{}

System::System(System const& s):
	n_(s.n_),
	N_(s.N_), 
	m_(s.m_),
	M_(s.M_), 
	bc_(s.bc_),
	links_(s.links_)
{}

System::~System(){}
/*}*/

void System::save(IOFiles& w) const {
	w("energy per site",E_);
	w("correlation on links",corr_);
	w("long range correlation",long_range_corr_);
}
