#include "MCParticle.hpp"

void MCParticle::init_Particle(double fx){
	Particle::init_Particle(fx);
	unsigned int s(history_.size());
	if(s){
		Rand<unsigned int> rnd(1,s);
		s = rnd.get();
		while(history_.target_next() && --s);
		x_ = history_.get().get_param();
	}
}

void MCParticle::print() const {
	Particle::print();
	std::cout<<"particle history ("<<history_.size()<<")"<<std::endl;
	history_.set_target();
	while( history_.target_next() ){
		std::cout<<history_.get_ptr()<<" ";
		history_.get().print();
		std::cout<<std::endl;
	}
}

bool MCParticle::update(std::shared_ptr<MCSim> const& new_elem){
	if(history_.find_sorted(new_elem,MCSim::cmp_for_merge)){ history_.set_target(); }
	else{ history_.add_after_target(new_elem); }

	/*\warning may not need to run select_new_best at each step*/
	if(Nupdate_ == update_now_){
		Nupdate_ = 0;
		return select_new_best();
	} else {
		Nupdate_++;
		double tmp(new_elem->get_S()->get_energy().get_x());
		if(tmp<fbx_){
			set_bx_via(new_elem->get_param());
			fbx_ = tmp;
			return true;
		}
		return false;
	}
}

void MCParticle::add_to_history(std::shared_ptr<MCSim> const& new_elem){
	history_.add_end(new_elem);
}

bool MCParticle::select_new_best(){
	double tmp;
	Vector<double> param;
	while(history_.target_next()){
		tmp = history_.get().get_S()->get_energy().get_x();
		if(tmp<fbx_){
			param = history_.get().get_param();
			fbx_ = tmp;
		}
	}
	if(param.ptr()){
		set_bx_via(param); 
		return true;
	} else {
		return false;
	}
}

void MCParticle::set_bx_via(Vector<double> const& param){
	bool found(false);
	for(unsigned int i(0);i<Nfreedom_;i++){
		for(unsigned int j(0);j<ps_[i].size();j++){
			if( my::are_equal(ps_[i](j),param(i)) ){
				bx_(i) = j;
				j = ps_[i].size();
				found = true;
			}
		}
	}
	assert(found);
	(void)(found);
}

Vector<double> MCParticle::get_param() const {
	Vector<double> param(Nfreedom_);
	for(unsigned int i(0);i<Nfreedom_;i++){
		param(i) = ps_[i](floor(x_(i)));
	}
	return param;
}
