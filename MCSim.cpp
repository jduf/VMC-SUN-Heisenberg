#include "MCSim.hpp"

/*{constructors, destructors*/
MCSim::MCSim(Vector<double> const& param):
	param_(param)
{}

MCSim::MCSim(IOFiles& r):
	param_(r)
{
	unsigned int ref_type_of_MCSystem(r.read<unsigned int>());
	switch(ref_type_of_MCSystem){
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
/*}*/

void MCSim::create_S(System const* const s){
	CreateSystem cs(s);
	cs.init(&param_,NULL);
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
		}
	}
	if(!is_created()){
		std::cerr<<"void MCSim::create_S(Container* C) : faulty parameters : "<<param_<<std::endl;
	}
}

void MCSim::copy_S(std::unique_ptr<MCSystem> const& S){
	S_ = S->clone();
}

void MCSim::run(unsigned int const& thermalization_steps, unsigned int const& tmax){
	if(is_created()){
		MonteCarlo mc(S_.get(),tmax);
		mc.thermalize(thermalization_steps);
		mc.run();
	} else {
		std::cerr<<"void MCSim::run(unsigned int const& thermalization_steps, unsigned int const& tmax) : faulty parameters : "<<param_<<std::endl;
	}
}

bool MCSim::check_conv(double const& convergence_criterion){
	S_->get_energy().complete_analysis(convergence_criterion);
	return S_->get_energy().get_conv();
}

void MCSim::complete_analysis(double const& convergence_criterion){
	S_->complete_analysis(convergence_criterion);
}

void MCSim::set_observable(bool all){
	S_->set_observable(all);
}

void MCSim::free_memory(){
	S_->free_memory();
}

/*{static methods*/
bool MCSim::sort_by_E(MCSim const& a, MCSim const& b){
	return a.get_S()->get_energy().get_x()<b.get_S()->get_energy().get_x();
}

unsigned int MCSim::sort_by_param_for_merge(MCSim const& list, MCSim const& new_elem){
	for(unsigned int i(0);i<list.param_.size();i++){
		if(list.param_(i) - new_elem.param_(i) > 0.0001){ return 0; }
		if(list.param_(i) - new_elem.param_(i) <-0.0001){ return 1; }
	}
	return 2;
}

void MCSim::merge(MCSim& list, MCSim& new_elem){
	list.get_S()->get_energy().merge(new_elem.get_S()->get_energy());
}
/*}*/

void MCSim::write(IOFiles& w) const {
	w<<param_<<S_->get_ref()(1);
	S_->write(w);
}

void MCSim::save(IOFiles& w) const {
	CreateSystem cs(S_.get());
	cs.init(&param_,NULL);
	if(cs.get_status()==2){
		cs.save_param(w);
		S_->save_input(w);
		S_->save_output(w);
	}
}
