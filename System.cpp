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
}

System::System(System const& S):
	n_(S.n_),
	N_(S.N_),
	m_(S.m_),
	M_(S.M_),
	bc_(S.bc_),
	links_(S.links_)
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

