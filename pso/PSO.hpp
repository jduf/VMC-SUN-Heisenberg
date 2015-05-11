#ifndef DEF_PSO
#define DEF_PSO

#include<memory>

#include"Rand.hpp"
#include"Vector.hpp"
#include"omp.h"

class Optimization{
	public:
		Optimization() = default;
		virtual ~Optimization() = default;
		/*{Forbidden*/
		Optimization(Optimization const&) = delete;
		Optimization(Optimization&&) = delete;
		Optimization& operator=(Optimization) = delete;
		/*}*/

		static void set_limit(unsigned int const& param, double const& min, double const& max);
		static void set(unsigned int const& Nfreedom, double const& cg, double const& cp);
		static bool within_limit(Vector<double> const& x);
		static unsigned int get_Nfreedom() { return Nfreedom_; }

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
		virtual ~Particle() = default;
		/*{Forbidden*/
		Particle(Particle const&) = delete;
		Particle(Particle&&) = delete;
		Particle& operator=(Particle) = delete;
		/*}*/

		virtual void init_Particle(double fx);
		virtual void move(Vector<double> const& bx_all);
		virtual void print() const;

		void set_best(Vector<double> const& x, double fx);
		double const& get_fbx() const { return fbx_; }
		Vector<double> const& get_x() const { return x_; }
		Vector<double> const& get_bx() const { return bx_; }
		Vector<double> const& get_v() const { return v_; }

	protected:
		double fbx_;		//!< value at the best position
		Vector<double> x_;	//!< position
		Vector<double> v_;	//!< velocity
		Vector<double> bx_;	//!< best position 

	private:
		Rand<double> rnd_;
};

template<typename Type>
class Swarm{
	public:
		Swarm(unsigned int const& Nparticles, unsigned int const& maxiter, unsigned int const& Nfreedom, double const& cg, double const& cp);
		~Swarm();
		/*{Forbidden*/
		Swarm() = delete;
		Swarm(Swarm const&) = delete;
		Swarm(Swarm&&) = delete;
		Swarm& operator=(Swarm const&) = delete;
		/*}*/

		void init_PSO(double const& fx);
		void run(double const& tol);
		void print() const;

		Vector<double> const& get_bx() const { return particle_[bparticle_]->get_bx(); }

	protected:
		unsigned int const Nparticles_;//!< numbre of particles
		std::vector<std::shared_ptr<Particle> > particle_;

	private:
		unsigned int const maxiter_;//!< maximum number of iteration
		unsigned int bparticle_;	//!< best particle
		bool* free_;         		//!< true if particle_ isn't running

		/*!This method must exist is the child class, it is the function that
		 * is minimized*/
		virtual bool evaluate(unsigned int const& p)=0;
		void next_step(unsigned int const& p);
		bool is_converged(double const& tol);
};

/*{Swarm*/
/*constructors and destructor*/
/*{*/
template<typename Type>
Swarm<Type>::Swarm(unsigned int const& Nparticles, unsigned int const& maxiter, unsigned int const& Nfreedom, double const& cg, double const& cp):
	Nparticles_(Nparticles),
	particle_(Nparticles),
	maxiter_(maxiter),
	bparticle_(0),
	free_(new bool[Nparticles])
{
	for(unsigned int i(0);i<Nparticles_;i++){ particle_[i] = std::make_shared<Type>();	}
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
void Swarm<Type>::init_PSO(double const& fx){
	for(unsigned int p(0);p<Nparticles_;p++){ free_[p] = true; }

#pragma omp parallel for schedule(dynamic,1)
	for(unsigned int p=0;p<Nparticles_;p++){
		particle_[p]->init_Particle(fx);
		evaluate(p);
	}
	/*as bparticle_=0, start at i=1*/
	for(unsigned int p(1);p<Nparticles_;p++){
		if(particle_[p]->get_fbx() < particle_[bparticle_]->get_fbx() ){ bparticle_ = p; }
	}
}

template<typename Type>
void Swarm<Type>::run(double const& tol){
	if(int(Nparticles_)<=omp_get_max_threads()){
		for(unsigned int i(0);i<maxiter_;i++){
#pragma omp parallel for
			for(unsigned int p=0;p<Nparticles_;p++){ next_step(p); }
			if(is_converged(tol)){ i=maxiter_; }
		}
	} else {
		bool keepon(true);
		unsigned int p(0);
		unsigned int local_p(0);
#pragma omp parallel for schedule(dynamic,1) private(local_p)
		for(unsigned int i=0; i<maxiter_*Nparticles_;i++){
#pragma omp critical
			{
				local_p=p;
				p = (p+1) % Nparticles_;
			}
			if(keepon){
				if(free_[local_p]){
					free_[local_p]=false;
					next_step(local_p);
					free_[local_p]=true;
					if(!((i+1)%Nparticles_)){ keepon = is_converged(tol); }
				} else { i--; }
			}
		}
	}
}

template<typename Type>
void Swarm<Type>::print() const {
	std::cout<<"Print each particle"<<std::endl;
	for(unsigned int i(0);i<Nparticles_;i++){
		particle_[i]->print();
	}
}
/*}*/

/*private methods*/
/*{*/
template<typename Type>
void Swarm<Type>::next_step(unsigned int const& p){
	particle_[p]->move(particle_[bparticle_]->get_x());
	if(evaluate(p)){
#pragma omp critical(update_best_pos_particle)
		{
			for(unsigned int p(0);p<Nparticles_;p++){
				if( particle_[p]->get_fbx() < particle_[bparticle_]->get_fbx() ){ bparticle_ = p; }
			}
		}
	}
}

template<typename Type>
bool Swarm<Type>::is_converged(double const& tol){
	Vector<double> dx(Nparticles_);
	Vector<double> bx(particle_[bparticle_]->get_bx());
	for(unsigned int j(0);j<Nparticles_;j++){
		dx(j) = ((particle_[j]->get_bx()-bx).norm_squared() + (particle_[j]->get_x()-bx).norm_squared())/Nparticles_;
	}
	return dx.mean() < tol;
}
/*}*/
/*}*/
#endif
