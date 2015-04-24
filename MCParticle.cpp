#include "MCParticle.hpp"

/*{MCParticle*/
void MCParticle::move(Vector<double> const& bx_all){
	Particle::move(bx_all);
	unsigned int n;
	double dx(0.01);
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

bool MCParticle::update(std::shared_ptr<MCSim> const& new_elem){
	if(history_.find_sorted(new_elem,MCSim::cmp_for_fuse)){ history_.set_free(); }
	else{ history_.add_after_free(new_elem); }

	double tmp;
	bool is_better(false);
	while(history_.go_to_next()){
		tmp = history_.get().get_S()->get_energy().get_x();
		if(tmp<fbx_){
			bx_ = history_.get().get_param();
			fbx_ = tmp;
			is_better = true;
		}
	}
	return is_better;
}

void MCParticle::print() const {
	Particle::print();
	std::cout<<"particle history "<<std::endl;
	history_.set_free();
	while( history_.go_to_next() ){
		std::cout<<history_.get_ptr()<<" ";
		history_.get().print();
	}
}
/*}*/
