#include "MCSim.hpp"

/*constructors, destructors*/
/*{*/
MCSim::MCSim(Vector<double> const& param):
	param_(param)
{}

MCSim::MCSim(IOFiles& r):
	param_(r)
{
	switch(r.read<unsigned int>()){
		case 0:
			{ MCS_ = std::unique_ptr<SystemBosonic<double> >(new SystemBosonic<double>(r)); } break;
		case 1:
			{ MCS_ = std::unique_ptr<SystemFermionic<double> >(new SystemFermionic<double>(r)); } break;
		case 2:
			{ MCS_ = std::unique_ptr<SystemFermionic<std::complex<double> > >(new SystemFermionic<std::complex<double> >(r)); } break;
	}
}
/*}*/

/*core methods*/
/*{*/
void MCSim::create_S(System const* const s){
	CreateSystem cs(s);
	cs.init(&param_,NULL);
	if(cs.get_status()==2){
		cs.create();
		if(cs.get_status()==1){
			if(cs.use_complex()){
				if(cs.is_bosonic()){
					MCS_.reset(new SystemBosonic<std::complex<double> >(*dynamic_cast<const Bosonic<std::complex<double> >*>(cs.get_GenericSystem())));
				} else {
					MCS_.reset(new SystemFermionic<std::complex<double> >(*dynamic_cast<const Fermionic<std::complex<double> >*>(cs.get_GenericSystem())));
				}
			} else {
				if(cs.is_bosonic()){
					MCS_.reset(new SystemBosonic<double>(*dynamic_cast<const Bosonic<double>*>(cs.get_GenericSystem())));
				} else {
					MCS_.reset(new SystemFermionic<double>(*dynamic_cast<const Fermionic<double>*>(cs.get_GenericSystem())));
				}
			}
		}
	}
	if(!is_created()){ std::cerr<<__PRETTY_FUNCTION__<<" : status_="<<cs.get_status()<<", faulty parameters="<<param_<<std::endl; }
}

void MCSim::copy_S(std::unique_ptr<MCSystem> const& S){
	MCS_ = S->clone();
}

void MCSim::run(unsigned int const& ts, unsigned int const& tmax){
	if(is_created()){
		MonteCarlo mc(MCS_.get(),tmax);
		mc.thermalize(ts);
		mc.run();
	} else {
		std::cerr<<__PRETTY_FUNCTION__<<" : faulty parameters : "<<param_<<std::endl;
	}
}
/*}*/

/*write in IOFiles methods*/
/*{*/
void MCSim::write(IOFiles& w) const {
	w<<param_<<MCS_->get_ref()(1);
	MCS_->write(w);
}

void MCSim::save(IOFiles& w) const {
	CreateSystem cs(MCS_.get());
	cs.init(&param_,NULL);
	if(cs.get_status()==2){ cs.save(w); }
	else { std::cerr<<__PRETTY_FUNCTION__<<" : status="<<cs.get_status()<<std::endl; }
}
/*}*/

/*static methods*/
/*{*/
bool MCSim::sort_by_E(MCSim const& a, MCSim const& b){
	return a.get_MCS()->get_energy().get_x()<b.get_MCS()->get_energy().get_x();
}

unsigned int MCSim::sort_for_merge(MCSim const& list, MCSim const& new_elem){
	return sort_by_param_for_merge(list.param_, new_elem.param_);
}

unsigned int MCSim::sort_by_param_for_merge(Vector<double> const& a, Vector<double> const& b){
	for(unsigned int i(0);i<a.size();i++){
		if(a(i) - b(i) > 0.0001){ return 0; }
		if(a(i) - b(i) <-0.0001){ return 1; }
	}
	return 2;
}
/*}*/
