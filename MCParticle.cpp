#include "MCParticle.hpp"

void MCParticle::move(Vector<double> const& bx_all){
	Particle::move(bx_all);
	/*!move to different parameter set can be achieved if v>1, therefore if
	 * the particle is static, it could be relaunched*/
	if(v_.norm_squared()<0.25){
		std::cerr<<__PRETTY_FUNCTION__<<" : static particle, reinitialization init_Particle(100)"<<std::endl;
		init_Particle(100);
	}
}

void MCParticle::print() const {
	Particle::print();
	std::cout<<"particle history ("<<history_.size()<<")"<<std::endl;
	history_.set_target();
	while( history_.target_next() ){
		std::cout<<history_.get_ptr()<<" "<<history_.get().get_param()<<" "<<history_.get().get_energy()<<std::endl;
	}
}

bool MCParticle::update(std::shared_ptr<MCSim> const& new_elem){
	if(!history_.find_in_sorted_list(new_elem,MCSim::sort_for_merge)){
		history_.add_after_target(new_elem);
	}

	/*\warning may not need to run select_new_best at each step*/
	if(Nupdate_ == update_now_){
		Nupdate_ = 0;
		return select_new_best();
	} else {
		Nupdate_++;
		double tmp(new_elem->get_energy().get_x());
		if(tmp<fbx_){
			set_bx_via(new_elem->get_param());
			fbx_ = tmp;
			return true;
		}
		return false;
	}
}

bool MCParticle::select_new_best(){
	double tmp;
	Vector<double> param;
	history_.set_target();
	while(history_.target_next()){
		tmp = history_.get().get_energy().get_x();
		if(tmp<fbx_){
			param = history_.get().get_param();
			fbx_ = tmp;
		}
	}
	if(param.ptr()){
		set_bx_via(param);
		return true;
	} else { return false; }
}

void MCParticle::set_bx_via(Vector<double> param){
	bool found;
	for(unsigned int i(0);i<dof_;i++){
		found = false;
		for(unsigned int j(0);j<ps_[i].size();j++){
			if( my::are_equal(ps_[i](j),param(i)) ){
				bx_(i) = j;
				j = ps_[i].size();
				found = true;
			}
		}
		if(!found){
			std::cerr<<__PRETTY_FUNCTION__<<" : can't find a match for dof "<<i<<" for param "<<param(i)<<" in "<<ps_[i]<<std::endl;
		}
		assert(found);
	}
}

Vector<double> MCParticle::get_param() const {
	Vector<double> param(dof_);
	for(unsigned int i(0);i<dof_;i++){
		if(x_(i)<min_(i) || x_(i)>=max_(i)){
			std::cerr<<__PRETTY_FUNCTION__<<" : "<<i<<" bug "<<x_<<" | "<<min_<<" | "<<max_<<std::endl;
			for(unsigned int j(0);j<dof_;j++){ std::cout<<ps_[j]<<std::endl; }
		}
		param(i) = ps_[i](floor(x_(i)));
	}
	return param;
}
