#include "System.hpp"

/*constructors and destructor*/
/*{*/
System::System(
		Vector<unsigned int> const& ref,
		unsigned int const& N, 
		unsigned int const& m,
		unsigned int const& n,
		int const& bc, 
		Vector<unsigned int> const& M,
		Vector<double> const& J):
	ref_(ref),
	N_(N), 
	m_(m),
	n_(n),
	bc_(bc),
	M_(M),
	J_(J),
	status_(4)
{
	if(M_.sum() != m_*n_ || m_>N_){ std::cerr<<"System::System(ref,N,n,m,bc,M,J) : Bad initialization"<<std::endl; } 
	else{status_--;}
	std::cout<<J_<<std::endl;
}

System::System(IOFiles& r):
	ref_(r),
	N_(r.read<unsigned int>()), 
	m_(r.read<unsigned int>()),
	n_(r.read<unsigned int>()),
	bc_(r.read<int>()),
	M_(r),
	J_(r),
	status_(r.read<unsigned int>()),
	links_(r),
	E_(r),
	corr_(r),
	lr_corr_(r)
{}
/*}*/

void System::set_observable(unsigned int const& what){
	E_.set(50,5,false);
	if(what>0){ corr_.set(links_.row(),50,5,false);}
	if(what>1){ lr_corr_.set(links_.row(),50,5,false); }
}

void System::set_binning(){ 
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
	w<<ref_<<N_<<m_<<n_<<bc_<<M_<<J_<<status_<<links_<<E_<<corr_<<lr_corr_;
}

void System::save(IOFiles& w) const {
	w.write("ref",ref_);
	w.write("N", N_);
	w.write("m", m_);
	w.write("n", n_);
	w.write("bc",bc_);
	w.write("M", M_);
	w.write("J", J_);
}

