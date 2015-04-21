#ifndef DEF_PSO
#define DEF_PSO

#include"Rand.hpp"
#include"Vector.hpp"
#include"omp.h"

#include<memory>

class Optimization{
	public:
		Optimization(){}
		virtual ~Optimization(){}

		static void set_limit(unsigned int const& param, double const& min, double const& max);
		static void set(unsigned int const& Nfreedom, double const& cg, double const& cp);

	protected:
		static unsigned int Nfreedom_; //!< number of degrees of freedom
		static Vector<double> min_;    //!< min_(c) minimum value of the cth coordinate
		static Vector<double> max_;    //!< max(c) minimum value of the cth coordinate
		static double cg_;	//! group influence
		static double cp_; //! personal influence
		static double chi_;//! constriction factor
};

class Particle: public Optimization{
	public:
		Particle():rnd_(0.0,1.0){}
		virtual ~Particle(){}

		void init(double fx);
		virtual void move(Vector<double> const& bx_all);

		double const& get_fbx() const { return fbx_; }
		Vector<double> const& get_x() const { return x_; }
		Vector<double> const& get_bx() const { return bx_; }
		void set_best(Vector<double> const& x, double fx){ 
			bx_ = x;
			fbx_ = fx; 
		}

		virtual void print(std::ostream& flux) const;

	protected:
		double fbx_;		//!< value at the best position
		Vector<double> x_;	//!< position
		Vector<double> v_;	//!< velocity
		Vector<double> bx_;	//!< best position 

	private:
		Particle(Particle const& p);

		Rand<double> rnd_;
};

std::ostream& operator<<(std::ostream& flux, Particle const& p);

template<typename Type>
class Swarm {
	public:
		Swarm(unsigned int const& Nparticles, unsigned int const& maxiter, unsigned int const& Nfreedom, double const& cg, double const& cp);
		~Swarm();

		void init(double const& fx);
		void run();
		void print(std::ostream& flux) const;

		Vector<double> const& get_bx() const { return p_[bparticle_]->get_bx(); }

	protected:
		std::vector<std::shared_ptr<Particle> > p_;

	private:
		unsigned int const Nparticles_;//!< numbre of particles
		unsigned int const maxiter_;	//!< maximum number of iteration
		unsigned int bparticle_;//!< best particle
		bool* free_;         //!< true if particle_ isn't running

		void next_step(unsigned int const& p);
		/*!This method must exist is the child class, it is the function that
		 * is minimized*/
		virtual bool is_better_x(unsigned int const& p)=0;
};

/*{Swarm*/
/*constructors and destructor*/
/*{*/
template<typename Type>
Swarm<Type>::Swarm(unsigned int const& Nparticles, unsigned int const& maxiter, unsigned int const& Nfreedom, double const& cg, double const& cp):
	p_(Nparticles),
	Nparticles_(Nparticles),
	maxiter_(maxiter),
	bparticle_(0),
	free_(new bool[Nparticles])
{
	for(unsigned int i(0);i<Nparticles_;i++){ p_[i] = std::make_shared<Type>();	}
	Optimization::set(Nfreedom,cg,cp);
}

template<typename Type>
Swarm<Type>::~Swarm(){
	if(free_){delete[] free_;}
}
/*}*/

/*public methods*/
/*{*/
template<typename Type>
void Swarm<Type>::init(double const& fx){
	for(unsigned int p(0);p<Nparticles_;p++){ free_[p] = true; }

#pragma omp parallel for schedule(dynamic,1)
	for(unsigned int p=0;p<Nparticles_;p++){
		p_[p]->init(fx);
		is_better_x(p);
	}
	/*as bparticle_=0, start at i=1*/
	for(unsigned int p(1);p<Nparticles_;p++){
		if(p_[p]->get_fbx() < p_[bparticle_]->get_fbx() ){ bparticle_ = p; }
	}
}

template<typename Type>
void Swarm<Type>::run(){
	if(int(Nparticles_)<=omp_get_max_threads()){
		std::cout<<"PSO::run all particles in parallel"<<std::endl;
		for(unsigned int iter(0);iter<maxiter_;iter++){
#pragma omp parallel for
			for(unsigned int p=0;p<Nparticles_;p++){ next_step(p); }
		}
	} else {
		std::cout<<"PSO::run only "+my::tostring(omp_get_max_threads())+" particles in parallel"<<std::endl;
		unsigned int p(0);
#pragma omp parallel for schedule(dynamic,1) firstprivate(p)
		for(unsigned int i=0; i<maxiter_*Nparticles_; i++){
			if(free_[p]){
				free_[p]=false;
				next_step(p);
				free_[p]=true;
			} else { i--; }
			p = (p+1) % Nparticles_;
		}
	}
}

template<typename Type>
void Swarm<Type>::print(std::ostream& flux) const {
	std::cout<<"Print each particle"<<std::endl;
	for(unsigned int i(0);i<Nparticles_;i++){
		p_[i]->print(flux);
	}
}
/*}*/

/*private methods*/
/*{*/
template<typename Type>
void Swarm<Type>::next_step(unsigned int const& p){
	p_[p]->move(p_[bparticle_]->get_x());
	if(is_better_x(p)){
#pragma omp critical(update_best_pos_particle)
		{
			for(unsigned int p(0);p<Nparticles_;p++){
				if( p_[p]->get_fbx() < p_[bparticle_]->get_fbx() ){ bparticle_ = p; }
			}
		}
	}
}
/*}*/
/*}*/
#endif
