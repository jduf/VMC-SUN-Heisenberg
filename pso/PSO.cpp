#include"PSO.hpp"

double Particle::cg_ = 2.05;
double Particle::cp_ = 2.05;
double Particle::chi_ = 0.792844;
unsigned int Particle::dof_ = 0;
Vector<double> Particle::min_ = Vector<double>();
Vector<double> Particle::max_ = Vector<double>();

void Particle::set(unsigned int const& dof, double const& cg, double const& cp){
	dof_ = dof;
	min_.set(dof,-1e10);
	max_.set(dof,1e10);
	chi_ = -2.0/(2.0-(cp_+cg_)-sqrt((cp_+cg_)*(cp_+cg_)-4*(cp_+cg_)));
	if(cg+cp>=4){
		cg_ = cg;
		cp_ = cp;
	} else {
		std::cerr<<__PRETTY_FUNCTION__<<" : cg+cp<4 => chi=nan. Keep definition cg=cp=2.05, chi=0.729844"<<std::endl;
	}
	chi_ = 2.0/(cp_+cg_-2.0+sqrt((cp_+cg_)*(cp_+cg_)-4*(cp_+cg_)));
}

void Particle::set_limit(unsigned int const& param, double const& min, double const& max){
	min_(param) = min;
	max_(param) = max;
}

Particle::Particle():
	rnd_(0.0,1.0),
	free_(true)
{}

void Particle::init_Particle(double fx){
	fbx_ = fx;
	x_.set(dof_);
	v_.set(dof_);

	for(unsigned int i(0);i<dof_;i++){
		x_(i) = rnd_.get()*(max_(i)-min_(i))+min_(i);
		v_(i) = rnd_.get()*(rnd_.get()>0.5?(max_(i)-x_(i)):(min_(i)-x_(i)));
	}
	bx_ = x_;
	move(bx_);
	bx_ = x_;
}

void Particle::move(Vector<double> const& bx_all){
	for(unsigned int i(0);i<dof_;i++){
		v_(i) = chi_*(v_(i) + cp_*rnd_.get()*(bx_(i)-x_(i)) + cg_*rnd_.get()*(bx_all(i)-x_(i)));
		if( x_(i)+v_(i) > max_(i)){ v_(i) = log(1.0+rnd_.get()*(exp(max_(i)-x_(i))-1.0)); }
		if( x_(i)+v_(i) < min_(i)){ v_(i) =-log(1.0+rnd_.get()*(exp(x_(i)-min_(i))-1.0)); }
		x_(i) = my::chop(x_(i)+v_(i));
		assert( x_(i) <= max_(i) );
		assert( x_(i) >= min_(i) );
	}
}

void Particle::print() const {
	std::cout<<"x="<<x_<<" bx="<<bx_<<" fbx="<<fbx_<<std::endl;
}
/*}*/
