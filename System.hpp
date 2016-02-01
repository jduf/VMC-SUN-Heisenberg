#ifndef DEF_SYSTEM
#define DEF_SYSTEM

#include "Parseur.hpp"
#include "Observable.hpp"

/*{*//*!Class that contains the parameters of the Hamiltonian and observables.
	   (SU(N), number of particles, bond strength,... and energy,
	   correlations,...)

	   Main class, everything relies on this class. It contains all the
	   important variables and results.

	   The status_ variables gives a way to stop a simulation if the
	   initialization failed and also indicates up to where it went :
	   + status_ = 5 : System is initialized (in System)
	   + status_ = 4 : System is constructed
	   + status_ = 3 : GenericSystem has a well defined unit cell size
	   + status_ = 2 : System*D has a good geometry
	   + status_ = 1 : Non degenerate eigenvalues (in System*D)
	   + status_ = 0 : Found an initial state (in SystemFermionic)
	   *//*}*/
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
		/*!Sets the observables to default (0) values and initilizes binning*/
		void set_obs(std::vector<Observable> const& obs, int const& nobs);
		/*!Remove all osbservable with index higher than from*/
		void clear_obs(unsigned int const& from);
		/*!Resets the observables*/
		void reset_obs();
		/*!Checks if the energy has converged to a stable value*/
		bool check_conv(double const& convergence_criterion);
		/*!Calls complete_analysis of the sampled datas*/
		void complete_analysis(double const& convergence_criterion);
		/*!Merge the binning of all observables of s into (*this)*/
		void merge(System* const s);
		/*!Deletes the binning for all observables*/
		void delete_binning();
		/*}*/

		/*{Saves/prints ref_, N_, m_, n_, bc_, M_, J_, status_ and obs_*/
		/*!Call child method to save the whole VMC without descritption*/
		virtual void write(IOFiles& w) const;
		/*!Saves only this class with description*/
		void save(IOFiles& w) const;
		/*!Prints the values of the nobs first observables*/
		void print(unsigned int nobs) const;
		/*}*/

		/*{Simple value return*/
		/*!Returns the reference to the type of wavefunction*/
		Vector<unsigned int> const& get_ref() const { return ref_; }
		/*!Returns the status of the system*/
		unsigned int const& get_status() const { return status_; }
		/*!Returns the energy*/
		Data<double> const& get_energy() const { return obs_[0][0]; }
		/*!Returns the energy*/
		unsigned int const& get_n() const { return n_; }
		/*!Returns the number of observables*/
		unsigned int nobs() const { return obs_.size(); }
		/*!Returns all observables*/
		std::vector<Observable> const& get_obs() const { return obs_; }
		/*}*/

	protected:
		/*!Default copy constructor*/
		System(System const&) = default;
		/*!Default constructor*/
		System():ref_(0),N_(0),m_(0),n_(0),bc_(0),M_(0),status_(5){ std::cout<<__PRETTY_FUNCTION__<<" : should never be called"<<std::endl; }

		Vector<unsigned int> const ref_;//!< type of system
		unsigned int const N_;			//!< number of colors
		unsigned int const m_;			//!< number of particles per site
		unsigned int const n_;			//!< number of sites
		int const bc_;					//!< boundary condition
		Vector<unsigned int> const M_;	//!< number of particles of each color

		/*!the following attributes will be set by GenericSystem and can't be
		 * defined within System*/
		Vector<double> J_;			//!< coupling strength
		unsigned int status_;		//!< status of the simulation
		std::vector<Observable> obs_;//!< observables (energy, bond energy, correlations, overlap...)

		/*!Saves ref_, N_, m_, n_, bc_, M_ and J_ in w (with description)*/
		void save_input(IOFiles& w) const;
		/*!Saves status_ and obs in w (with description)*/
		void save_output(IOFiles& w) const;

	private:
		/*!Complete the required parameters that are not yet in Parseur*/
		Vector<unsigned int> complete_system_info(Parseur& P);
};
#endif
