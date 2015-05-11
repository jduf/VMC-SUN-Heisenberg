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

unsigned int MCSim::cmp_for_fuse(MCSim const& list, MCSim const& new_elem) {
	for(unsigned int i(0);i<list.param_.size();i++){
		if(list.param_(i) > new_elem.param_(i)){ return 0; }
		if(list.param_(i) < new_elem.param_(i)){ return 1; }
	}
	return 2;
}

void MCSim::fuse(MCSim& list_elem, MCSim& new_elem) { 
	list_elem.get_S().get()->get_energy().merge(new_elem.get_S().get()->get_energy());
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
}

void MCSim::copy_S(std::unique_ptr<MCSystem> const& S){
	S_ = S->clone();
}

void MCSim::write(IOFiles& w) const {
	w<<ref_<<param_;
	S_->write(w);
}

void MCSim::save(Container* C) const {
	CreateSystem cs(C,&param_);
	cs.init();
	if(cs.get_status()==2){
		cs.create();
		if(cs.get_status()==1){
			Linux command;
			command("/bin/mkdir -p " + cs.get_path());
			IOFiles file_results(cs.get_path() + cs.get_filename()+".jdbin",true);
			cs.init_output_file(file_results);
			cs.save();
			RST rst;
			rst.title("Results",'-');
			file_results.add_header()->add(rst.get());
			file_results.write("energy per site",S_->get_energy());
			file_results.write("correlation on links",S_->get_corr());
			file_results.write("long range correlation",S_->get_lr_corr());
		}
	}
}
