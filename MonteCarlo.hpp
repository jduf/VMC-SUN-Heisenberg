#ifndef DEF_MONTECARLO
#define DEF_MONTECARLO

#include "SystemBosonic.hpp"
#include "SystemFermionic.hpp"

/*{! Class MonteCarlo
 *
 * Implement the Monte-Carlo algorithm. This alorithm let a System evolves
 * according to a Markov process. To work properly, System has to contain at
 * least these methods : 
 *
 * - System::ratio() : compute the probability to accept the next configuration
 * - System::update() : update the old cufiguration to the new one
 * - System::measure() : measure an observable according to the current configuration 
 *
 * The MonteCarlo class contains a pointer to System. According to the system,
 * the constructor creates a SystemFermionic or SystemBosonic. As the
 * constructors of these classes need a pointer to a CreateSystem, the
 * constructor of MonteCarlo must provide a pointer on the same class.
 *
 * Each time that a class is instanciated, a random number generator is created
 * according to the thread on which the code is running. The same thread number
 * is then transmitted to create the System.
 *
 * Once the System is created, it is thermalized.
 *
 * The MonteCarlo::run() method lunches the Monte-Carlo simulation.
 }*/
template <typename Type>
class MonteCarlo{
	public:
		/*!Constructor*/
		MonteCarlo(CreateSystem* CS, unsigned int tmax, unsigned int type); 
		/*!Simple destructor*/
		~MonteCarlo();

		/*!Run the Monte-Carlo algorithm*/
		void run();
		/*!Get the pointer on the system*/
		System<Type>* get_system() const { return S_;}

		void check();

	private:
		/*!Forbids the copy constructor*/
		MonteCarlo(MonteCarlo const& mc); 
		/*!Forbids the assignment operator*/
		MonteCarlo const& operator=(MonteCarlo const& mc);

		/*!Find the next configuration and measure it*/
		void next_step();
		//{Private method that gives a shutoff condition
		/*!Stops the simulation when
		 *
		 * - convergence is reached
		 * - time limit is up
		 * - if the energy diverges
		 *
		 * If the E_ diverges, the simulation is restarted
		 */
		//}
		bool keepon(double const& tol);

		unsigned int const tmax_;//!< Time limit in second, by default 5min
		System<Type>* S_; 		//!< Pointer to a Fermionic or Bosonic System 
		Rand* rnd_;				//!< Pointer to a random number generator
		Time time_; 			//!< To stop the simulation after time_limit seconds
};

/*constructors and destructor*/
/*{*/
template<typename Type>
MonteCarlo<Type>::MonteCarlo(CreateSystem* CS, unsigned int tmax, unsigned int type):
	tmax_(tmax),
	S_(NULL),
	rnd_(NULL)
{
	unsigned int thread(omp_get_thread_num());
	rnd_ = new Rand(1e4,thread);
	if(CS->is_bosonic()){S_ = new SystemBosonic<Type>(CS,thread,type);}
	else{S_ = new SystemFermionic<Type>(CS,thread,type);}
	if(S_->ready()){
		double ratio(0.0);
		for(unsigned int i(0);i<1e5;i++){
			S_->swap();
			ratio = norm_squared(S_->ratio());
			if( ratio > rnd_->get() ){ S_->update(); }
		}
		S_->measure_new_step();
	}
}

template<typename Type>
MonteCarlo<Type>::~MonteCarlo(){
	delete S_;
	delete rnd_;
}
/*}*/

/*public methods*/
/*{*/
template<typename Type>
void MonteCarlo<Type>::run(){
	if(S_->ready()){
		do{next_step();}
		while(keepon(5e-5));
	}
	S_->complete_analysis(5e-5);
}
/*}*/

/*private methods*/
/*{*/
template<typename Type>
void MonteCarlo<Type>::next_step(){
	S_->swap();
	if( norm_squared(S_->ratio()) > rnd_->get() ){
		S_->update();
		S_->measure_new_step();
	}
	S_->add_sample();
}

template<typename Type>
bool MonteCarlo<Type>::keepon(double const& tol){
	if(time_.limit_reached(tmax_)){ return false; }
	if(S_->is_converged(tol)){}
	//if(std::abs(S_->get_energy())>1e2){ 
		//std::cerr<<"Simulation diverges => is restarted"<<std::endl;
		//S_->set();
	//}
	return true;
}
/*}*/

template<typename Type>
void MonteCarlo<Type>::check(){
	unsigned int i(0);
	if(S_->ready()){/*passed the first two steps*/
		do{ i++; next_step();}
		while(keepon(1e-5));
	}
}
#endif
