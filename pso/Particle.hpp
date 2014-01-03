#ifndef DEF_PARTICLE
#define DEF_PARTICLE

#include"Rand.hpp"
#include"Functions.hpp"
#include<vector>
#include<iostream>

template<typename Type>
class Particle {
	public:
		Particle(unsigned int Nparam, double cg, double cp, Type const& Class, double (Type::*f)(std::vector<double>));
		~Particle();

		void init(unsigned int thread);
		void move(std::vector<double> gb);
		void evaluate();		
		void set_limit(unsigned int param, double min, double max);
		double get_fb(){return fb_;}
		std::vector<double> get_b(){return b_;}

		std::vector<double> path_fx;

		void result();
		void print();
		void path();

	private:
		std::vector<double> x_; //!< position
		std::vector<double> b_; //!< personal best position
		std::vector<double> v_; //!< velocity
		std::vector<double> min_; //!< velocity
		std::vector<double> max_; //!< velocity
		double cg_;
		double cp_;
		double chi_;
		double fb_; //!< value at the best position
		std::vector<std::vector<double> > path_;
		Rand* rnd_;
		Type Class_;
		double (Type::*f_)(std::vector<double>);
};

template<typename Type>
Particle<Type>::Particle(unsigned int Nparam, double cg, double cp, Type const& Class, double (Type::*f)(std::vector<double>)):
	x_(Nparam),
	b_(Nparam),
	v_(Nparam),
	min_(Nparam),
	max_(Nparam),
	cg_(cg),
	cp_(cp),
	chi_(-2.0/(2.0-(cp+cg)-sqrt((cp+cg)*(cp+cg)-4*(cp+cg)))),
	fb_(2),
	rnd_(NULL),
	Class_(Class),
	f_(f)
{
	if(cg_+cp_<4){
		std::cerr<<"Particle::Particle(Nparam,cg,cp,*f) : cg+cp<4 => chi=nan. Redefinition cg=cp=2.05"<<std::endl;
		cg_ = 2.05;
		cp_ = 2.05;
	}
}

template<typename Type>
Particle<Type>::~Particle(){
	if(rnd_){ delete rnd_; }
}

template<typename Type>
void Particle<Type>::init(unsigned int thread){
	rnd_ = new Rand(1e6,thread);
	for(unsigned int i(0);i<x_.size();i++){
		//min_[i] = -M_PI/2.0; 
		//max_[i] = 2.0*M_PI/3.0; 
		min_[i] = -1.5; 
		max_[i] = 1.5; 
		x_[i] = rnd_->get()*(max_[i]-min_[i]);
		v_[i] = rnd_->get()*(max_[i]-min_[i]);
		if(rnd_->get()>0.5){
			x_[i] = -x_[i];
		}
		if(rnd_->get()>0.5){
			v_[i] = -v_[i];
		}
	}
	fb_= (Class_.*f_)(x_);
	b_ = x_;
}

template<typename Type>
void Particle<Type>::evaluate(){
	double fx((Class_.*f_)(x_));
	//std::cout<<fx<<std::endl;
	if( fx < fb_){
		fb_ = fx;
		b_ = x_;
	}
	path_fx.push_back(fx);
}

template<typename Type>
void Particle<Type>::move(std::vector<double> gb){
	for(unsigned int i(0);i<x_.size();i++){
		v_[i] = chi_*(v_[i] + cp_*rnd_->get()*(b_[i]-x_[i]) + cg_*rnd_->get()*(gb[i]-x_[i]));
		//if(v_[i] > (max_[i]-min_[i])/2.0){v_[i] = (max_[i]-min_[i])/4.0;}
		//if(v_[i] < (min_[i]-max_[i])/2.0){v_[i] = (min_[i]-max_[i])/4.0;}
		if( x_[i]+v_[i] > max_[i]){ 
			v_[i] = log(1.0+rnd_->get()*(exp(max_[i]-x_[i])-1.0));
		}
		if( x_[i]+v_[i] < min_[i]){ 
			v_[i] = -log(1.0+rnd_->get()*(exp(x_[i]-min_[i])-1.0));
		}
		x_[i] += v_[i]; 
	}
	path_.push_back(x_);
}

template<typename Type>
void Particle<Type>::set_limit(unsigned int param, double min, double max){
	min_[param] = min;
	max_[param] = max;
}

template<typename Type>
void Particle<Type>::print(){
	for(unsigned int i(0);i<x_.size();i++){
		std::cout<<x_[i]<<" ";
	}
	std::cout<<std::endl;
}

template<typename Type>
void Particle<Type>::path(){
	for(unsigned int i(0);i<path_.size();i++){
		std::cout<<i<<" ";
		for(unsigned int j(0);j<x_.size();j++){
			std::cout<<path_[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

template<typename Type>
void Particle<Type>::result(){
	std::cerr<<"best particle ";
	for(unsigned int j(0);j<x_.size();j++){
		std::cerr<<b_[j]<<" ";
	}
	std::cerr<<" -> "<<fb_<<std::endl;
}
#endif
