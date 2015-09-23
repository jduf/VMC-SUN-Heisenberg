#ifndef DEF_MONTECARLO
#define DEF_MONTECARLO

#include "SystemBosonic.hpp"
#include "SystemFermionic.hpp"

/*{!Class MonteCarlo
 *
 * Implement the Monte-Carlo algorithm. This alorithm lets a System evolve
 * according to a Markov process. To work properly, MCSystem has to contain at
 * least these methods : 
 *
 * - Type System::ratio() : compute the probability to accept the next
 *   configuration
 * - void System::update() : update the old cufiguration to the new one
 * - void System::measure() : measure an observable according to the current
 *   configuration 
 *
 * The MonteCarlo class contains a pointer to MCSystem. The MCSystem is in
 * reality either a SystemFermionic or a SystemBosonic. 
 *
 * Once the System is created, it is thermalized.
 *
 * - void MonteCarlo::run() : lunch the Monte-Carlo simulation.
 * - void MonteCarlo::complete_analysis(double tol) : make sure that the
 *   computed quantities are fully computed
 }*/
class MonteCarlo{
	public:
		/*!Constructor*/
		MonteCarlo(MCSystem* S, unsigned int const& tmax);
		/*!Default destructor*/
		~MonteCarlo() = default;
		/*{Forbidden*/
		MonteCarlo() = delete; 
		MonteCarlo(MonteCarlo const&) = delete; 
		MonteCarlo(MonteCarlo&&) = delete; 
		MonteCarlo& operator=(MonteCarlo) = delete; 
		/*}*/

		/*!Thermalize the Monte-Carlo algorithm*/
		void thermalize(unsigned int const& thermalization_steps);
		/*!Run the Monte-Carlo algorithm*/
		void run();

	private:
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
		double ratio_;			 //!< Ratio between current and next step
		MCSystem* S_;			 //!< Pointer to a Fermionic or Bosonic System 
		Time time_; 			 //!< To stop the simulation after time_limit seconds
		Rand<double> rnd_;		 //!< Pointer to a random number generator
};
#endif
