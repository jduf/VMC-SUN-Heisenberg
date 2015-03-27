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
	min_.set(Nfreedom,-2.5);
	max_.set(Nfreedom,2.5);
	chi_ = -2.0/(2.0-(cp_+cg_)-sqrt((cp_+cg_)*(cp_+cg_)-4*(cp_+cg_)));
	if(cg+cp>=4){
		cg_ = cg;
		cp_ = cp;
	} else {
		std::cerr<<"Optimization::set(unsigned int Nfreedom, double cg, double cp) : cg+cp<4 => chi=nan. Keep definition cg=cp=2.05, chi=0.729844"<<std::endl;
	}
	chi_ = -2.0/(2.0-(cp_+cg_)-sqrt((cp_+cg_)*(cp_+cg_)-4*(cp_+cg_)));
}

void Optimization::set_limit(unsigned int const& param, double const& min, double const& max){
	min_(param) = min;
	max_(param) = max;
}
/*}*/

/*{Particle*/
/*core of the class*/
/*{*/
void Particle::move(Vector<double> const& bx_all){
	for(unsigned int i(0);i<Nfreedom_;i++){
		v_(i) = chi_*(v_(i) + cp_*rnd_.get()*(bx_(i)-x_(i)) + cg_*rnd_.get()*(bx_all(i)-x_(i)));
		if( x_(i)+v_(i) > max_(i)){ v_(i) = log(1.0+rnd_.get()*(exp(max_(i)-x_(i))-1.0)); }
		if( x_(i)+v_(i) < min_(i)){ v_(i) =-log(1.0+rnd_.get()*(exp(x_(i)-min_(i))-1.0)); }
		x_(i) += v_(i); 
	}
}

void Particle::init(double fx){
	fbx_ = fx;
	x_.set(Nfreedom_,1.0);
	v_.set(Nfreedom_,2.0);

	for(unsigned int i(0);i<Nfreedom_;i++){
		x_(i) = rnd_.get()*(max_(i)-min_(i))+min_(i);
		v_(i) = rnd_.get()*(max_(i)-min_(i))+min_(i);
	}

	bx_ = x_;
	move(x_);
	bx_ = x_;
}
/*}*/
/*}*/

