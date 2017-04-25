#include "VMCAnalyse.hpp"

VMCAnalyse::VMCAnalyse(IOFiles& in):
	VMCMinimization(in,true,"ANA")
{}

VMCAnalyse::VMCAnalyse(List<IOFiles>& in):
	VMCMinimization(in.first(),true,"MERGE")
{
	std::string info_tmp(m_->info_.get());
	m_->info_.set();
	m_->info_.title("Merge "+my::tostring(in.size())+" different simulations",'+');

	m_->info_.text(info_tmp);
	update_info(m_->samples_);

	in.set_target();
	in.target_next();
	while(in.target_next()){
		std::cout<<"will merge "<<in.get_ptr()->get_filename()<<std::endl;
		VMCAnalyse merging(in.get());
		m_->info_.text(merging.m_->info_.get());
		update_info(merging.m_->samples_);

		unsigned int iter(0);
		merging.m_->samples_.set_target();
		while(merging.m_->samples_.target_next()){
			my::display_progress(iter++,merging.m_->samples_.size(),"merging");
			if( m_->samples_.find_in_sorted_list(merging.m_->samples_.get_ptr(),MCSim::sort_for_merge) ){ m_->samples_.get().merge(merging.m_->samples_.get_ptr()); }
			else { m_->samples_.add_after_target(merging.m_->samples_.get_ptr()); }
		}
	}
}

void VMCAnalyse::update_info(List<MCSim> const& merging_samples){
	double E(0);
	List<MCSim>::Node* best(NULL);
	merging_samples.set_target();
	while(merging_samples.target_next()){
		if(merging_samples.get().get_energy().get_x()<E){
			E = merging_samples.get().get_energy().get_x();
			best = merging_samples.get_target();
		}
	}
	std::string p("(");
	for(unsigned int i(0);i<best->get()->get_param().size()-1;i++){ p += my::tostring(best->get()->get_param()(i))+","; }
	p += my::tostring(best->get()->get_param().back())+")";
	m_->info_.def("Best parameter","p="+p);
	m_->info_.def("Best energy","E="+my::tostring(best->get()->get_energy().get_x())+", dEoE="+my::tostring(best->get()->get_energy().get_dx()));
	m_->info_.np();
}
