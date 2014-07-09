#include "System.hpp"

/*constructors and destructor*/
/*{*/
System::System(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc):
	ref_(ref),
	N_(N), 
	m_(m),
	n_(n),
	M_(M),
	bc_(bc),
	status_(3)
{
	if(M_.sum() != m_*n_ || m_>N_){ std::cerr<<"System::System(N,n,m,M,bc,ref) : Bad initialization"<<std::endl; } 
	else{status_--;}
}

System::System(System const& s):
	ref_(s.ref_),
	N_(s.N_), 
	m_(s.m_),
	n_(s.n_),
	M_(s.M_),
	bc_(s.bc_),
	status_(s.status_),
	E_(s.E_),
	corr_(s.corr_),
	links_(s.links_)
{}
/*}*/

void System::set(){ 
	E_.set(); 
	corr_.set(); 
	long_range_corr_.set(); 
}
