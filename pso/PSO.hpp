#ifndef DEF_PSO
#define DEF_PSO

#include <memory>
#include "Rand.hpp"
#include "Vector.hpp"
#include "omp.h"

class Particle{
	public:
		Particle();
		virtual ~Particle() = default;
		/*{Forbidden*/
		Particle(Particle const&) = delete;
		Particle(Particle&&) = delete;
		Particle& operator=(Particle const&) = delete;
		/*}*/

		virtual void init_Particle(double fx);
		virtual void move(Vector<double> const& bx_all);
		virtual void print() const;

		bool is_free() const { return free_; }
		void toggle_free(){ free_ = !free_; }
		double const& get_fbx() const { return fbx_; }
		Vector<double> const& get_x() const { return x_; }
		Vector<double> const& get_bx() const { return bx_; }
		Vector<double> const& get_v() const { return v_; }

		static void set_limit(unsigned int const& param, double const& min, double const& max);
		static void set(unsigned int const& dof, double const& cg, double const& cp);

		void define_current_as_best(double const& fbx){ fbx_ = fbx; bx_ = x_; }
		
	protected:
		double fbx_;		//!< value at the best position
		Vector<double> x_;	//!< position
		Vector<double> v_;	//!< velocity
		Vector<double> bx_;	//!< best position
		static Vector<double> min_;//!< min_(c) minimum value of the cth coordinate
		static Vector<double> max_;//!< max(c) minimum value of the cth coordinate
		static unsigned int dof_;  //!< number of degrees of freedom

	private:
		Rand<double> rnd_;
		bool free_;

		static double cg_;	//!< group influence
		static double cp_;  //!< personal influence
		static double chi_; //!< constriction factor
};

template<typename Type>
class Swarm{
	public:
		Swarm(unsigned int const& Nparticles, unsigned int const& maxiter, unsigned int const& dof, double const& cg, double const& cp);
		~Swarm() = default;
		/*{Forbidden*/
		Swarm() = delete;
		Swarm(Swarm const&) = delete;
		Swarm(Swarm&&) = delete;
		Swarm& operator=(Swarm const&) = delete;
		/*}*/

		void init_PSO(double const& fx);
		void minimize();
		void print() const;

		Vector<double> const& get_bx() const { return particle_[bparticle_]->get_bx(); }

	protected:
		unsigned int const Nparticles_;//!< numbre of particles
		unsigned int const maxiter_;   //!< maximum number of iteration
		std::vector<std::shared_ptr<Particle> > particle_;

	private:
		unsigned int bparticle_;//!< best particle

		/*!This method must exist in the child class, it is the function that
		 * is minimized*/
		virtual bool evaluate(unsigned int const& p)=0;
		void next_step(unsigned int const& p);
};

/*{Swarm*/
/*constructors and destructor*/
/*{*/
template<typename Type>
Swarm<Type>::Swarm(unsigned int const& Nparticles, unsigned int const& maxiter, unsigned int const& dof, double const& cg, double const& cp):
	Nparticles_(Nparticles),
	maxiter_(maxiter),
	particle_(Nparticles),
	bparticle_(0)
{
	for(unsigned int i(0);i<Nparticles_;i++){ particle_[i] = std::make_shared<Type>(); }
	Particle::set(dof,cg,cp);
}
/*}*/

/*public methods*/
/*{*/
template<typename Type>
void Swarm<Type>::init_PSO(double const& fx){
#pragma omp parallel for
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
void Swarm<Type>::minimize(){
	if(int(Nparticles_)<=omp_get_max_threads()){
		for(unsigned int i(0);i<maxiter_;i++){
#pragma omp parallel for
			for(unsigned int p=0;p<Nparticles_;p++){ next_step(p); }
		}
	} else {
		unsigned int p(0);
		unsigned int local_p(0);
#pragma omp parallel for schedule(dynamic) private(local_p)
		for(unsigned int i=0;i<maxiter_*Nparticles_;i++){
#pragma omp critical(Swarm__minimize__local)
			{
				local_p=p;
				p = (p+1) % Nparticles_;
				if(particle_[local_p]->is_free()){ particle_[local_p]->toggle_free(); }
				else { local_p = Nparticles_; }
			}
			if(local_p<Nparticles_){
				next_step(local_p);/*!will call the function to minimize*/
				particle_[local_p]->toggle_free();
			} else { i--; }
		}
	}
}

template<typename Type>
void Swarm<Type>::print() const {
	std::cout<<"Print each particle"<<std::endl;
	for(unsigned int i(0);i<Nparticles_;i++){ particle_[i]->print(); }
}
/*}*/

/*private methods*/
/*{*/
template<typename Type>
void Swarm<Type>::next_step(unsigned int const& p){
	particle_[p]->move(particle_[bparticle_]->get_x());
	if(evaluate(p)){
#pragma omp critical(Swarm__next_step__local)
		{
			for(unsigned int p(0);p<Nparticles_;p++){
				if( particle_[p]->get_fbx() < particle_[bparticle_]->get_fbx() ){ bparticle_ = p; }
			}
		}
	}
}
/*}*/
/*}*/
#endif
