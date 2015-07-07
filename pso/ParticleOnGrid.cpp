#include "ParticleOnGrid.hpp"

bool Measure::func(Measure const& a, Measure const& b) { 
	unsigned int i(0);
	while(i<a.x_.size()){
		if(a.x_(i) > b.x_(i)){ return false; }
		if(a.x_(i) < b.x_(i)){ return true; }
		if(a.x_(i)== b.x_(i)){ i++; }
	}
	return false;
}

std::ostream& operator<<(std::ostream& flux, Measure const& m){
	m.print(flux);
	return flux;
}

ParticleOnGrid::~ParticleOnGrid(){
#pragma omp critical
	{
		std::cout<<particle_history_<<std::endl;
	}
}

void ParticleOnGrid::move(Vector<double> const& bx_all){
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

void ParticleOnGrid::update(double fx){
	bx_ = x_;
	fbx_ = fx;
}

void ParticleOnGrid::add_history(std::shared_ptr<Measure> const& m){
	particle_history_.add_sort(m,Measure::func);
}
