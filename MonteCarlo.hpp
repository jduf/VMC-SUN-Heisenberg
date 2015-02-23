#ifndef DEF_MONTECARLO
#define DEF_MONTECARLO

#include "SystemBosonic.hpp"
#include "SystemFermionic.hpp"
#include <omp.h>

/*{! Class MonteCarlo
 *
 * Implement the Monte-Carlo algorithm. This alorithm lets a System evolve
 * according to a Markov process. To work properly, MCSystem has to contain at
 * least these methods : 
 *
 * - Type System::ratio() : compute the probability to accept the next configuration
 * - void System::update() : update the old cufiguration to the new one
 * - void System::measure() : measure an observable according to the current configuration 
 *
 * The MonteCarlo class contains a pointer to MCSystem. The MCSystem is in
 * reality either a SystemFermionic or a SystemBosonic. 
 *
 * Each time that a class is instanciated, a random number generator is created
 * according to the seeds another random nunber generator provides. The same
 * thread number is then transmitted to create the System.
 *
 * Once the System is created, it is thermalized.
 *
 * The MonteCarlo::run() method lunches the Monte-Carlo simulation.
 }*/
template <typename Type>
class MonteCarlo{
	public:
		/*!Constructor*/
		MonteCarlo(MCSystem<Type>* S, unsigned int const& tmax);
		/*!Simple destructor*/
		~MonteCarlo(){}

		/*!Thermalize the Monte-Carlo algorithm*/
		void thermalize(unsigned int const& N);
		/*!Run the Monte-Carlo algorithm*/
		void run();
		void complete_analysis(double tol){ S_->complete_analysis(tol); }
		void delete_binning(){ S_->delete_binning();}

	private:
		/*!Forbids copy*/
		MonteCarlo(MonteCarlo const& mc); 
		/*!Forbids assignment*/
		MonteCarlo const& operator=(MonteCarlo const& mc);

		/*!Find the next configuration and measure it*/
		void next_step();
		/*{Description*/
		/*!Private method that gives a shutoff condition 
		 * Stops the simulation when
		 * - convergence is reached
		 * - time limit is up
		 * - kill=true
		 *
		 * If the E_ diverges, the simulation is restarted */
		/*}*/
		bool keepon();

		unsigned int const tmax_;//!< Time limit in second, by default 5min
		MCSystem<Type>* S_;		//!< Pointer to a Fermionic or Bosonic System 
		Time time_; 			//!< To stop the simulation after time_limit seconds
		Rand<double> rnd_;		//!< Pointer to a random number generator
};

/*constructors and destructor*/
/*{*/
template<typename Type>
MonteCarlo<Type>::MonteCarlo(MCSystem<Type>* S, unsigned int const& tmax):
	tmax_(tmax),
	S_(S),
	rnd_(0.0,1.0)
{}
/*}*/

/*public methods*/
/*{*/
template<typename Type>
void MonteCarlo<Type>::thermalize(unsigned int const& N){
	if(S_->get_status()==0){
		for(unsigned int i(0);i<N;i++){
			S_->swap();
			if( norm_squared(S_->ratio()) > rnd_.get() ){ S_->update(); }
		}
		S_->measure_new_step();
	}
}

template<typename Type>
void MonteCarlo<Type>::run(){
	if(S_->get_status()==0){
		do{next_step();}
		while(keepon());
	}
}
/*}*/

/*private methods*/
/*{*/
template<typename Type>
void MonteCarlo<Type>::next_step(){
	S_->swap();
	if( norm_squared(S_->ratio()) > rnd_.get() ){
		S_->update();
		S_->measure_new_step();
	}
	S_->add_sample();
}

template<typename Type>
bool MonteCarlo<Type>::keepon(){
	if(time_.limit_reached(tmax_)){ return false; }
	if(time_.progress(tmax_/21)){
		if(!omp_get_thread_num()){
			S_->get_energy().compute_convergence(1e-5);
			std::cerr<<"E="<<S_->get_energy().get_x()<<" ("<<S_->get_energy().get_dx()<<") after "<<100.0*time_.elapsed()/tmax_<<"%"<<std::endl;
			//S_->set();
		}
	}
	if(std::abs(S_->get_energy().get_x())>1e2){ 
		std::cerr<<"Simulation diverges (E="<<S_->get_energy().get_x()<<") => is restarted"<<std::endl;
		S_->set();
	}
	return true;
}
/*}*/
#endif
