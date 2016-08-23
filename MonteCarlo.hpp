#ifndef DEF_MONTECARLO
#define DEF_MONTECARLO

#include "SystemBosonic.hpp"
#include "SystemFermionic.hpp"
#include "SystemBiFermionic.hpp"
#include <omp.h>

/*{*//*!Class that implement the Monte-Carlo importance sampling algorithm on
	   configuration stored in MCSystem.

	   This alorithm lets a System evolve according to a Markov process. To
	   work properly, MCSystem has to contain at least these methods :

	   + Type System::ratio() : compute the probability to accept the next
	   configuration
	   + void System::update() : update the old cufiguration to the new one
	   + void System::measure() : measure an observable according to the
	   current configuration

	   The MonteCarlo class contains a pointer to MCSystem. The MCSystem is in
	   reality a child of MCSystem (e.g. SystemFermionic).

	   Once the System is created, it should be measured as follows :

	   + void MonteCarlo::thermalize(unsigned int const& ts) : thermalization
	   phase to reach configuration of relevant importance.
	   + void MonteCarlo::run() : lunch the Monte-Carlo simulation.
	   *//*}*/
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
		void thermalize(unsigned int const& ts);
		/*!Run the Monte-Carlo algorithm*/
		void run();
		/*!Run the Monte-Carlo algorithm*/
		void run(unsigned int const& maxiter);

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
		unsigned int iter_;		 //!< internal counter (keepon() check each iter=1e5)
		double ratio_;			 //!< Ratio between current and next step
		MCSystem* S_;			 //!< Pointer to a Fermionic or Bosonic System
		Time time_; 			 //!< To stop the simulation after time_limit seconds
		Rand<double> rnd_;		 //!< Pointer to a random number generator
};
#endif
