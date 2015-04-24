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
	status_(4)
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
	lr_corr_(s.lr_corr_),
	links_(s.links_)
{}

System::System(IOFiles& r):
	ref_(r),
	N_(r.read<unsigned int>()), 
	m_(r.read<unsigned int>()),
	n_(r.read<unsigned int>()),
	M_(r),
	bc_(r.read<int>()),
	status_(r.read<unsigned int>()),
	E_(r),
	corr_(r),
	lr_corr_(r),
	links_(r)
{}
/*}*/

void System::set(){ 
	E_.set(); 
	corr_.set(); 
	lr_corr_.set(); 
}

void System::delete_binning(){ 
	E_.delete_binning();
	corr_.delete_binning();
	lr_corr_.delete_binning();
}

void System::write(IOFiles& w) const {
	w<<ref_<<N_<<m_<<n_<<M_<<bc_<<status_<<E_<<corr_<<lr_corr_<<links_;
}
