#ifndef DEF_SWARM
#define DEF_SWARM

#include"Particle.hpp"
#include"Time.hpp"
#include<omp.h>

template<typename Type>
class Swarm {
	public:
		Swarm(unsigned int Nparticle, unsigned int Nparam, double cg, double cp, Type const& Class, double (Type::*f)(std::vector<double>));
		
		void set_limit(unsigned int param, double min, double max);
		void launch(unsigned int i);
		bool converged() const;
		void next_step();
		void run();

		void result();
		void print();
		void path();

	private:
		std::vector<Particle<Type> > particle_;
		std::vector<bool> free_;//!< true if particle_[i] isn't running
		std::vector<double> b_; //!< global best position
		double fb_; 			//!< value at the global best position
		unsigned int maxiter_;
		unsigned int iter_;
		bool converged_;

		/*what follows will be deleted*/
		Time t_;
		std::vector<std::vector<time_t> > history_;
};

template<typename Type>
Swarm<Type>::Swarm(unsigned int Nparticle, unsigned int Nparam, double cg, double cp, Type const& Class, double (Type::*f)(std::vector<double>)):
	particle_(Nparticle,Particle<Type>(Nparam,cg,cp,Class,f)),
	free_(Nparticle,true),
	b_(Nparam),
	fb_(0),
	maxiter_(20),
	iter_(0),
	converged_(false),
	history_(Nparticle,std::vector<time_t>(2*maxiter_))
{
#pragma omp parallel for
	for(unsigned int i=0;i<particle_.size();i++){
		particle_[i].init(omp_get_thread_num());
		particle_[i].evaluate();
	}
	fb_= particle_[0].get_fb();
	b_ = particle_[0].get_b();
	history_[0][iter_] = t_.elapsed();
	for(unsigned int i(1);i<particle_.size();i++){
		if(particle_[i].get_fb() < fb_ ){
			fb_= particle_[i].get_fb();
			b_ = particle_[i].get_b();
		}
		history_[i][iter_] = t_.elapsed();
	}
}

template<typename Type>
void Swarm<Type>::set_limit(unsigned int param, double min, double max){
	for(unsigned int i(0);i<particle_.size();i++){
		particle_[i].set_limit(param,min,max);
	}
}

template<typename Type>
void Swarm<Type>::launch(unsigned int i){
	particle_[i].move(b_);
	particle_[i].evaluate();
	if(particle_[i].get_fb() < fb_ ){
		fb_ = particle_[i].get_fb();
		b_ = particle_[i].get_b();
	}
}

template<typename Type>
bool Swarm<Type>::converged() const {
	return converged_;
}

template<typename Type>
void Swarm<Type>::next_step(){
#pragma omp parallel for
	for(unsigned int i=0;i<particle_.size();i++){
		launch(i);
		history_[i][iter_] = t_.elapsed();
	}
}

template<typename Type>
void Swarm<Type>::run(){
	std::cerr<<"particle_.size() must be bigger than omp_get_num_threads"<<std::endl;
	unsigned int i(0);
	unsigned int ip(0);
#pragma omp parallel for schedule(dynamic,1) firstprivate(ip)
	for(unsigned int iter=0; iter<maxiter_; iter++){
#pragma omp critical
		{
			i++;
			ip = (i-1) % particle_.size();
			std::cout<<iter<<" "<<i<<" "<<ip;
			if(!free_[ip]){ 
				iter -= 1; 
			} else {
				std::cout<<" launched";
			}
			std::cout<<std::endl;
		}
		if(free_[ip]){
			free_[ip]=false;
			if(ip<2){sleep(2);}
			else{ sleep(5);}
			//launch(ip);
			free_[ip]=true;
		}
	}
}

template<typename Type>
void Swarm<Type>::result(){
	for(unsigned int i(0);i<b_.size();i++){
		std::cerr<<b_[i]<<" ";
	}
	std::cerr<<" -> "<<fb_<<std::endl;
	//for(unsigned int i(0);i<particle_.size();i++){
	//particle_[i].result();
	//}
}

template<typename Type>
void Swarm<Type>::print(){
	for(unsigned int i(0);i<particle_.size();i++){
		particle_[i].print();
	}
}

template<typename Type>
void Swarm<Type>::path(){
	for(unsigned int i(0);i<particle_.size();i++){
		for(unsigned int j(0);j<maxiter_;j++){
			std::cout<<history_[i][j]<<" "<<particle_[i].path_fx[j]<<std::endl;
		}
		std::cout<<std::endl;
	}
}
#endif
