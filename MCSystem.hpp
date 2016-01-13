#ifndef DEF_MCSYSTEM
#define DEF_MCSYSTEM

#include "System.hpp"
#include "Rand.hpp"
#include <memory>

/*{*//*!Abstract class that is used by MonteCarlo to sample the System 
	   using the Variational Monte-Carlo algorithm.

	   This class makes the connection between MonteCarlo and more specialized
	   System like SystemFermionic, SystemBiFermionic and SystemBosonic. This
	   is the reason why it is an abstract class with many (pure) virtual
	   methods.*//*}*/
class MCSystem: public virtual System{
	public:
		/*!Constructor*/
		MCSystem(System const& S);
		/*!Constructor that reads from file*/
		MCSystem(IOFiles& r);
		/*!Default destructor*/
		virtual ~MCSystem() = default;
		/*{Forbidden*/
		MCSystem() = delete;
		MCSystem(MCSystem&&) = delete;
		MCSystem& operator=(MCSystem const&) = delete;
		/*}*/

		/*!Exchanges two particles of different colors on random sites*/
		virtual void swap();
		/*!Exchanges the p0's particle of site s0 with the p1's of site s1*/
		virtual void swap(unsigned int const& s0, unsigned int const& s1, unsigned int const& p0, unsigned int const& p1);
		/*!Pure virtual method that computes the ratio between two states*/
		virtual double ratio()=0;
		/*!Updates only s_*/
		virtual void update();
		/*!Measures the system for the new step*/
		virtual void measure_new_step();
		/*!Adds the sample to the statistic*/
		void add_sample();

		/*!Returns a copy of the instance of the relevant child class*/
		virtual std::unique_ptr<MCSystem> clone() const = 0;
		/*!Saves some RAM*/
		virtual void free_memory() = 0;
		/*!Writes the curent state of the system, the color configuration, the
		 * observables and everything relevent to the simulation*/
		virtual void write(IOFiles& w) const;

	protected:
		/*!Allows copy if called from child class*/
		MCSystem(MCSystem const& mcsim);

		unsigned int new_c_[2];//!< colors of the exchanged sites
		unsigned int new_s_[2];//!< sites that are exchanged
		unsigned int new_p_[2];//!< particle on site that are exchanged

		Matrix<unsigned int> s_;  //!< s_(site,particle)=color

	private:
		/*!Checks only if the new state has not the same color on one site*/
		bool is_new_state_forbidden();

		Rand<unsigned int> n_rnd_;//!< generator of random numbers
		Rand<unsigned int> m_rnd_;//!< generator of random numbers

		unsigned int bei_;
};
#endif
