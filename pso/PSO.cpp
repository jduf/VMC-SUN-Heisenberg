#include"PSO.hpp"

/*{Optimization*/
double Optimization::cg_ = 2.05;
double Optimization::cp_ = 2.05;
double Optimization::chi_ = 0.792844;
unsigned int Optimization::Nfreedom_ = 0;
Vector<double> Optimization::min_ = Vector<double>();
Vector<double> Optimization::max_ = Vector<double>();

void Optimization::set(unsigned int const& Nfreedom, double const& cg, double const& cp){
	Nfreedom_ = Nfreedom;
	min_.set(Nfreedom,-1e10);
	max_.set(Nfreedom,1e10);
	chi_ = -2.0/(2.0-(cp_+cg_)-sqrt((cp_+cg_)*(cp_+cg_)-4*(cp_+cg_)));
	if(cg+cp>=4){
		cg_ = cg;
		cp_ = cp;
	} else {
		std::cerr<<"Optimization::set(unsigned int Nfreedom, double cg, double cp) : cg+cp<4 => chi=nan. Keep definition cg=cp=2.05, chi=0.729844"<<std::endl;
	}
	chi_ = 2.0/(cp_+cg_-2.0+sqrt((cp_+cg_)*(cp_+cg_)-4*(cp_+cg_)));
}

void Optimization::set_limit(unsigned int const& param, double const& min, double const& max){
	min_(param) = min;
	max_(param) = max;
}

bool Optimization::within_limit(Vector<double> const& x){
	for(unsigned int i(0);i<Nfreedom_;i++){
		if(x(i)<min_(i) || x(i)>max_(i)) { return false; }
	}
	return true;
}
/*}*/

/*{Particle*/
void Particle::init_Particle(double fx){
	fbx_ = fx;
	x_.set(Nfreedom_);
	v_.set(Nfreedom_);

	for(unsigned int i(0);i<Nfreedom_;i++){
		x_(i) = rnd_.get()*(max_(i)-min_(i))+min_(i);
		v_(i) = rnd_.get()*(max_(i)-min_(i))+min_(i);
	}
	bx_ = x_;
	move(bx_);
	bx_ = x_;
}

void Particle::move(Vector<double> const& bx_all){
	double r;
	for(unsigned int i(0);i<Nfreedom_;i++){
		v_(i) = chi_*(v_(i) + cp_*rnd_.get()*(bx_(i)-x_(i)) + cg_*rnd_.get()*(bx_all(i)-x_(i)));
		r = rnd_.get();
		if( x_(i)+v_(i) > max_(i)){ v_(i) = log(1.0+r*(exp(max_(i)-x_(i))-1.0)); }
		if( x_(i)+v_(i) < min_(i)){ v_(i) =-log(1.0+r*(exp(x_(i)-min_(i))-1.0)); }
#pragma omp critical
		{
			if( x_(i)+v_(i) > max_(i)){ std::cerr<<"there upper limit "<<x_(i)<<" "<<v_(i)<<" "<<max_(i)<<" "<<r<<std::endl; }
			if( x_(i)+v_(i) < min_(i)){ std::cerr<<"there lower limit "<<x_(i)<<" "<<v_(i)<<" "<<min_(i)<<" "<<r<<std::endl; }
		}
		x_(i) += v_(i); 
	}
}

void Particle::print() const {
	std::cout<<"x="<<x_<<" bx="<<bx_<<" fbx="<<fbx_<<std::endl;
}

void Particle::set_best(Vector<double> const& x, double fx){ 
	bx_ = x;
	fbx_ = fx; 
}
/*}*/
