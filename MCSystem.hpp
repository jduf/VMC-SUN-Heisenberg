#ifndef DEF_MCSYSTEM
#define DEF_MCSYSTEM

#include "System.hpp"
#include "Rand.hpp"
#include <memory>

/*!Abstract class that is used by MonteCarlo.hpp to sample the system.*/
class MCSystem: public virtual System{
	public:
		/*!Constructor*/
		MCSystem(System const& S);
		/*!Destructor*/
		virtual ~MCSystem(){}

		/*!Exchanges two particles of different colors on random sites*/
		virtual void swap();
		/*!Exchanges the p0's particle of site s0 with the p1's of site s1*/
		virtual void swap(unsigned int const& s0, unsigned int const& s1, unsigned int const& p0, unsigned int const& p1);
		/*!Pure virtual method that computes the ratio between two states*/
		virtual double ratio(bool const& squared)=0;
		/*!Updates only s_*/
		virtual void update();

		/*!Sample the system for the new step*/
		void measure_new_step();
		/*!Add the sample to the statistic*/
		void add_sample();
		/*!Calls complete_analysis of the sampled datas*/
		void complete_analysis(double const& tol);

		virtual std::unique_ptr<MCSystem> clone() const = 0;
		
	protected:
		unsigned int new_c_[2];//!< colors of the exchanged sites
		unsigned int new_s_[2];//!< sites that are exchanged
		unsigned int new_p_[2];//!< sites that are exchanged

		Matrix<unsigned int> s_;  //!< s(site,particle)=color
		Rand<unsigned int> n_rnd_;//!< generator of random numbers 
		Rand<unsigned int> m_rnd_;//!< generator of random numbers 

		/*!Forbid copy*/
		MCSystem(MCSystem const& mcsim);

	private:
		/*!Forbids default*/
		MCSystem();
		/*!Forbid assigment*/
		MCSystem& operator=(MCSystem const& mcsim);

		/*!Check only if the new state has not the same color on one site*/
		bool is_new_state_forbidden();
};
#endif
