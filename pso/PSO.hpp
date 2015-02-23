#ifndef DEF_PSO
#define DEF_PSO

#include"Rand.hpp"
#include"Vector.hpp"
#include"omp.h"

class PSO {
	public:
		PSO(unsigned int Nparticles, unsigned int Nfreedom, double cg, double cp, unsigned int maxiter);
		virtual ~PSO();

		void PSO_init();
		void PSO_set_limit(unsigned int param, double min, double max);
		void PSO_run();
		void PSO_print();
		void PSO_save(std::string filename);
		void PSO_load(std::string filename);

private:
		unsigned int Nparticles_;//!< numbre of particles
		unsigned int bparticle_;//!< best particle
		unsigned int Nfreedom_; //!< number of degrees of freedom
		unsigned int maxiter_;	//!< maximum number of iteration
		Vector<double> min_;    //!< min_(c) minimum value of the cth coordinate
		Vector<double> max_;    //!< max(c) minimum value of the cth coordinate
		Vector<double>* pbx_;	//!< best position 
		Vector<double>* pv_;	//!< velocity
		Vector<double>* px_;	//!< position
		bool* free_;            //!< true if particle_[i] isn't running
		double* pfbx_;			//!< value at the best position
		double cg_;	//! group influence
		double cp_; //! personal influence
		double chi_;//! constriction factor
		Rand<double> rnd_;

		void next_step(unsigned int i);
		void move(unsigned int i);
		void move_on_grid(unsigned int i);
		void evaluate(unsigned int i);		

		/*!This method must exist is the child class, it is the function that
		 * is minimized*/
		virtual double f(Vector<double> const& x)=0;
};
#endif
