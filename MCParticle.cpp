#include "MCParticle.hpp"

void MCParticle::move(Vector<double> const& bx_all){
	Particle::move(bx_all);
	unsigned int n;
	double dx(0.1);
	for(unsigned int j(0);j<Nfreedom_;j++){
		n=0;
		if(std::abs(x_(j))<dx/2){ n=1; }
		if(std::abs(x_(j)-min_(j))<dx/2){ n=2; }
		if(std::abs(x_(j)-max_(j))<dx/2){ n=3; }
		switch(n){
			case 0:{ x_(j) = std::round(x_(j)/dx)*dx; }break;
			case 1:{ x_(j) = 0; }break;
			case 2:{ x_(j) = min_(j); }break;
			case 3:{ x_(j) = max_(j); }break;
		}
	}
}

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
	if(history_.find_sorted(new_elem,MCSim::cmp_for_fuse)){ history_.set_target(); }
	else{ history_.add_after_target(new_elem); }

	/*\warning may not need to run select_new_best at each step*/
	if(Nupdate_ == update_now_){
		Nupdate_ = 0;
		return select_new_best();
	} else {
		Nupdate_++;
		double tmp(new_elem->get_S()->get_energy().get_x());
		if(tmp<fbx_){
			bx_ = new_elem.get()->get_param();
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
	bool new_best(false);
	while(history_.target_next()){
		tmp = history_.get().get_S()->get_energy().get_x();
		if(tmp<fbx_){
			bx_ = history_.get().get_param();
			fbx_ = tmp;
			new_best = true;
		}
	}
	return new_best;
}
