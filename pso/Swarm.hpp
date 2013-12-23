#ifndef DEF_SWARM
#define DEF_SWARM

#include"Particle.hpp"
#include<omp.h>

unsigned int MAXITER(10000);

class Swarm {
	public:
		Swarm(unsigned int Nparticle,unsigned int Nparam, double cg, double cp, double (*f)(std::vector<double>));
		
		void set_limit(unsigned int param, double min, double max);
		void next_step();
		void result();
		bool check();

		void print();
		void path();

	private:
		Swarm(const& Swarm);

		std::vector<Particle> particle_;
		std::vector<double> b_; //!< global best position
		double fb_; //!< value at the global best position
		unsigned int iter_;
};

Swarm::Swarm(unsigned int Nparticle,unsigned int Nparam, double cg, double cp, double (*f)(std::vector<double>)):
	particle_(Nparticle,Particle(Nparam,cg,cp,f)),
	b_(Nparam),
	fb_(2),
	iter_(0)
{
#pragma omp parallel for
	for(unsigned int i=0;i<particle_.size();i++){
		particle_[i].init(omp_get_thread_num());
		particle_[i].evaluate();
		if(particle_[i].get_fb() < fb_ ){
			fb_ = particle_[i].get_fb();
			b_ = particle_[i].get_b();
		}
	}
}

void Swarm::set_limit(unsigned int param, double min, double max){
	for(unsigned int i(0);i<particle_.size();i++){
		particle_[i].set_limit(param,min,max);
	}
}

void Swarm::next_step(){
#pragma omp parallel for
	for(unsigned int i=0;i<particle_.size();i++){
		particle_[i].move(b_);
		particle_[i].evaluate();
		if(particle_[i].get_fb() < fb_ ){
			fb_ = particle_[i].get_fb();
			b_ = particle_[i].get_b();
		}
	}
}

bool Swarm::check(){
	iter_++;
	if(iter_>MAXITER){return false; }
	else{ return true;}
}

void Swarm::result(){
	for(unsigned int j(0);j<b_.size();j++){
		std::cerr<<b_[j]<<" ";
	}
	std::cerr<<" -> "<<fb_<<std::endl;
	//for(unsigned int i(0);i<particle_.size();i++){
	//particle_[i].result();
	//}
}

void Swarm::print(){
	for(unsigned int i(0);i<particle_.size();i++){
		particle_[i].print();
	}
}

void Swarm::path(){
	for(unsigned int i(0);i<particle_.size();i++){
		particle_[i].path();
	}
	std::cout<<std::endl;
}
#endif
