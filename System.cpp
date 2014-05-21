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

System::System():
	n_(0),
	N_(0), 
	m_(0),
	M_(0), 
	bc_(0)
{
	std::cout<<"ok default System"<<std::endl;
}

System::~System(){}

void System::save(IOFiles& w) const {
	std::cout<<E_<<std::endl;
	w("energy per site",E_);
	w("correlation on links",corr_);
	w("long range correlation",long_range_corr_);
}
/*}*/

