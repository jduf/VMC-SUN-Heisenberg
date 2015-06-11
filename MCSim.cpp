#include "MCSim.hpp"

MCSim::MCSim(Vector<double> const& param):
	param_(param)
{}

MCSim::MCSim(IOFiles& r):
	ref_(r),
	param_(r)
{
	switch(ref_(1)){
		case 0:
			{
				S_ = std::unique_ptr<SystemBosonic<double> >(new SystemBosonic<double>(r));
			} break;
		case 1:
			{
				S_ = std::unique_ptr<SystemFermionic<double> >(new SystemFermionic<double>(r));
			} break;
		case 2:
			{
				S_ = std::unique_ptr<SystemFermionic<std::complex<double> > >(new SystemFermionic<std::complex<double> >(r));
			} break;
	}
}

unsigned int MCSim::cmp_for_fuse(MCSim const& list, MCSim const& new_elem){
	for(unsigned int i(0);i<list.param_.size();i++){
		if(list.param_(i) > new_elem.param_(i)){ return 0; }
		if(list.param_(i) < new_elem.param_(i)){ return 1; }
	}
	return 2;
}

void MCSim::fuse(MCSim& list, MCSim& new_elem) { 
	list.get_S()->get_energy().merge(new_elem.get_S()->get_energy());
}

void MCSim::create_S(Container* C, Vector<double> const* param){
	CreateSystem cs(C,param);
	cs.init();
	if(cs.get_status()==2){
		cs.create();
		if(cs.get_status()==1){
			if( cs.use_complex()){
				if(cs.is_bosonic()){
					S_.reset(new SystemBosonic<std::complex<double> >(*dynamic_cast<const Bosonic<std::complex<double> >*>(cs.get_system())));
				} else {
					S_.reset(new SystemFermionic<std::complex<double> >(*dynamic_cast<const Fermionic<std::complex<double> >*>(cs.get_system())));
				}
			} else {
				if(cs.is_bosonic()){
					S_.reset(new SystemBosonic<double>(*dynamic_cast<const Bosonic<double>*>(cs.get_system())));
				} else {
					S_.reset(new SystemFermionic<double>(*dynamic_cast<const Fermionic<double>*>(cs.get_system())));
				}
			}
			ref_ = cs.get_ref();
		}
	}
	if(!is_created()){
		std::cerr<<"void MCSim::create_S(Container* C) : faulty parameters : "<<param_<<std::endl;
	}
}

void MCSim::copy_S(std::unique_ptr<MCSystem> const& S){
	S_ = S->clone();
}

void MCSim::write(IOFiles& w) const {
	w<<ref_<<param_;
	S_->write(w);
}
