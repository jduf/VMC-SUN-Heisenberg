#include "System.hpp"

/*constructors and destructor and initialization*/
/*{*/
System::System(unsigned int N, unsigned int n, unsigned int m, int bc):
	n_(n),
	N_(N), 
	m_(m),
	M_((m_*n_)/N_), 
	bc_(bc)
{
	std::cout<<"ok Nnmbc System"<<std::endl; 
}

System::System(System const& s):
	n_(s.n_),
	N_(s.N_), 
	m_(s.m_),
	M_(s.M_), 
	bc_(s.bc_),
	links_(s.links_)
{
	std::cout<<"ok copy System"<<std::endl; 
}

System::~System(){}

void System::save(IOFiles& w) const {
	std::cout<<E_<<std::endl;
	w("energy per site",E_);
	w("correlation on links",corr_);
	w("long range correlation",long_range_corr_);
}
/*}*/

