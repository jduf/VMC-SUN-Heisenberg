#include "MCParticle.hpp"

void MCParticle::move(Vector<double> const& bx_all){
	Vector<double> tmp(x_);
	Particle::move(bx_all);
	unsigned int n;
	double dx(0.1);
	for(unsigned int i(0);i<Nfreedom_;i++){
		n=0;
		if(std::abs(x_(i))<dx/2){ n=1; }
		if(std::abs(x_(i)-min_(i))<dx/2){ n=2; }
		if(std::abs(x_(i)-max_(i))<dx/2){ n=3; }
		switch(n){
			case 0:{ x_(i) = std::round(x_(i)/dx)*dx; }break;
			case 1:{ x_(i) = 0; }break;
			case 2:{ x_(i) = min_(i)+dx; }break;
			case 3:{ x_(i) = max_(i)-dx; }break;
		}
	}
#pragma omp critical
	if(my::are_equal(tmp,x_)){
		Rand<int> rnd(0,2);
		std::cerr<<"force move"<<std::endl;
		double dv;
		for(unsigned int i(0);i<Nfreedom_;i++){
			dv = dx*rnd.get();
			v_(i) += dv;
			x_(i) += dv;
		}
	} else {
		std::cerr<<tmp<<" "<<x_<<std::endl;
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
		double tmp(new_elem.get()->get_S()->get_energy().get_x());
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
