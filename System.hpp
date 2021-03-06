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
	   + status_ = 1 : No degenerate eigenvalues (in System*D) or multiple identical observable set
	   + status_ = 0 : Found an initial state (in SystemFermionic)

	   The ref_ variable references the wavefunction, cluster shape...:
	   + ref_(0) : 1=chain
	               2=ladder
				   3=triangle
				   4=square
				   5=kagome
				   6=honeycomb
	   + ref_(1) : 1=real wave function
	               2=complex wave function
	   + ref_(2) : kind of wavefunction
	   + ref_(3) : geometry of the cluster
	   + ref_(4) : 0=do not create the cluster
	               1=create the cluster
				   2=create the cluster and the energy observable
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
		/*!Sets an other cluster geometry*/
		bool try_other_geometry(Vector<unsigned int> const& ref) const;
		/*!Sets coupling terms*/
		void set_J(System const* const s){ J_ = s->J_; }

		/*!Sets an observable*/
		void set_obs(Observable const& obs){ obs_.push_back(obs); }
		/*!Removes osbservables from 'from'*/
		void clear_obs(unsigned int const& from){ obs_.erase(obs_.begin()+from,obs_.end()); }
		/*!Reinitilizes the binnings (clear any previous measure)*/
		void reset_obs(){ for(auto& o:obs_){ o.reset(); } }

		/*!Checks if the energy has converged*/
		bool check_conv();
		/*!Calls complete_analysis on the observables*/
		void complete_analysis();
		/*!Merges the binning of all observables of s into (*this)*/
		void merge(System* const s);
		/*!Deletes the binning for all observables*/
		void delete_binning();
		/*}*/

		/*{Saves/prints ref_, N_, m_, n_, bc_, M_, J_, status_ and obs_*/
		/*!Call child method to save MCSystem's children without description*/
		virtual void write(IOFiles& w) const;
		/*!Saves this class and GenericSystem's children with description*/
		void save(IOFiles& w) const;
		/*!Prints the results of the simulation in std::cout*/
		void print(bool const& all) const;
		/*}*/

		/*{Simple value return*/
		/*!Returns the reference to the type of wavefunction*/
		Vector<unsigned int> const& get_ref() const { return ref_; }
		/*!Returns the status of the system*/
		unsigned int const& get_status() const { return status_; }
		/*!Returns the energy*/
		Data<double> const& get_energy() const { return obs_[0][0]; }
		/*!Gets relative incertitude of the energy*/
		double get_dEoE(){ return obs_[0][0].get_dx()/obs_[0][0].get_x(); }
		/*!Returns the energy*/
		unsigned int const& get_n() const { return n_; }
		/*!Returns the coupling strength*/
		Vector<double> const& get_J() const { return J_; }
		/*!Returns the number of observables*/
		unsigned int nobs() const { return obs_.size(); }
		/*!Returns all observables*/
		std::vector<Observable> const& get_obs() const { return obs_; }
		/*}*/

		void create_cluster(bool const& create){ ref_(4) = create; }

	protected:
		/*!Default copy constructor*/
		System(System const&) = default;
		/*!Default constructor*/
		System():ref_(0),N_(0),m_(0),n_(0),bc_(0),M_(0),status_(5){ std::cout<<__PRETTY_FUNCTION__<<" : should never be called"<<std::endl; }

		mutable Vector<unsigned int> ref_;//!< type of system
		unsigned int const N_;			//!< number of colors
		unsigned int const m_;			//!< number of particles per site
		unsigned int const n_;			//!< number of sites
		int const bc_;					//!< boundary condition
		Vector<unsigned int> const M_;	//!< number of particles of each color

		/*!The following attributes will be set by GenericSystem and can't be
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
