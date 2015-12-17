#ifndef DEF_SYSTEM
#define DEF_SYSTEM

#include "Parseur.hpp"
#include "Observable.hpp"

/*!{Class that contains the information on the state*/
/*
 * status_ = 5 : System is initialized (in System)
 * status_ = 4 : System is allowed (in System)
 * status_ = 3 : Boundary condition is allowed
 * status_ = 2 : Degenerated wavefunction
 * status_ = 1 : System*D is allowed (in System*D)
 * status_ = 0 : Found an initial state (in SystemFermionic)
 */
/*}*/
class System{
	public:
		/*!Constructor*/
		System(Parseur& P);
		/*!Constructor that reads from file*/
		System(IOFiles& r);
		/*!Default destructor*/
		virtual ~System() = default;
		/*{Forbidden*/
		System(System&&) = delete;
		System& operator=(System) = delete;
		/*}*/

		/*{Handles class attributes*/
		void set_observables(std::vector<Observable> const& obs, int const& nobs);
		/*!Sets the observables to default (0) values and initilizes binning*/
		void clear_observables(int const& nobs);
		/*!Checks if the energy has converged to a stable value*/
		bool check_conv(double const& convergence_criterion);
		/*!Calls complete_analysis of the sampled datas*/
		void complete_analysis(double const& convergence_criterion);
		/*!Merge the binning of all observables of s into (*this)*/
		void merge(System* const s);
		/*!Deletes the binning for all observables*/
		void delete_binning();
		/*}*/

		/*{Write in IOFiles methods*/
		virtual void write(IOFiles& w) const;
		/*{Description*/
		/*!Saves ref_, N_, m_, n_, bc_, M_ and J_ in w. The file will contain a
		 * brief description for each variable and eventually a header. As the
		 * method is virtual, a call on this method will call first
		 * child::save() const if it exists*/
		/*}*/
		virtual void save_input(IOFiles& w) const;
		void save_output(IOFiles& w) const;
		/*}*/

		/*{Simple value return*/
		/*!Returns the reference to the type of wavefunction*/
		Vector<unsigned int> const& get_ref() const { return ref_; }
		/*!Returns the status of the system*/
		unsigned int const& get_status() const { return status_; }
		/*!Returns the energy*/
		Data<double> const& get_energy() const { return E_; }
		/*!Returns all observables*/
		std::vector<Observable> const& get_obs() const { return obs_; }
		/*}*/

	protected:
		/*!Default copy constructor*/
		System(System const&) = default;
		/*!Default constructor*/
		System():ref_(0),N_(0),m_(0),n_(0),bc_(0),M_(0),status_(5){ std::cout<<__PRETTY_FUNCTION__<<" : should never be called"<<std::endl; }

		Vector<unsigned int> set_ref(Parseur& P);

		Vector<unsigned int> const ref_;//!< type of system
		unsigned int const N_;			//!< number of colors
		unsigned int const m_;			//!< number of particles per site
		unsigned int const n_;			//!< number of sites
		int const bc_;					//!< boundary condition
		Vector<unsigned int> const M_;	//!< number of particles of each color

		/*!the following attributes will be set by GenericSystem and can't be
		 * defined within System*/
		Vector<double> J_;				//!< coupling strength
		unsigned int status_;			//!< status of the simulation
		Data<double> E_; 				//!< energy of the system
		std::vector<Observable> obs_;	//!< all other observables (bond energy, correlations...)
};
#endif
