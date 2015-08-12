#ifndef DEF_SYSTEM
#define DEF_SYSTEM

#include "Sampling.hpp"
#include "Parseur.hpp"

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
		System(System const&) = delete;
		System(System&&) = delete;
		System& operator=(System) = delete;
		/*}*/

		std::vector<std::string> names() const;

		/*!Returns the reference to the type of wavefunction*/
		Vector<unsigned int> const& get_ref() const { return ref_; }
		/*!Returns the energy*/
		Data<double> const& get_energy() const { return E_; }
		/*!Returns a ref to E_ (needed for merging)*/
		Data<double>& get_energy(){ return E_; }
		/*!Returns the correlations*/
		DataSet<double> const& get_corr() const { return corr_; }
		/*!Returns the long range correlations*/
		DataSet<double> const& get_lr_corr() const { return lr_corr_; }
		/*!Returns the status of the system*/
		unsigned int const& get_status() const { return status_; }

		void set_bonds(System const* const s);
		/*!Sets the observables to default (0) values and initilizes binning*/
		void set_binning();
		/*!Sets the binning for E_(>=0), corr_(>=1), lr_corr(>=2)*/
		void set_observable(unsigned int const& what);
		/*!Deletes the binning for all observables*/
		void delete_binning();

		virtual void write(IOFiles& w) const;
		/*{Description*/
		/*!Saves ref_, N_, m_, n_, M_ bc_ and J_ in w. The file will contain a
		 * brief description for each variable and eventually a header. As the
		 * method is virtual, a call on this method will call first
		 * child::save() const if it exists*/
		/*}*/
		virtual void save_input(IOFiles& w) const;
		void save_output(IOFiles& w) const;

	protected:
		/*!Almost Copy constructor*/
		System(System const& s, unsigned int const& status);
		/*!Default constructor*/
		System():ref_(0),N_(0),m_(0),n_(0),bc_(0){ std::cout<<"System::System() : should never be called"<<std::endl; }

		Vector<unsigned int> const ref_;//!< type of system 
		unsigned int const N_;			//!< number of colors
		unsigned int const m_;			//!< number of particles per site
		unsigned int const n_;			//!< number of sites
		int const bc_;					//!< boundary condition
		Vector<unsigned int> const M_;	//!< number of particles of each color 
		Vector<double> J_;				//!< bond energy

		unsigned int status_		= 5;//!< status of the simulation
		Matrix<unsigned int> links_;	//!< bond <i,j>

		Data<double> E_; 				//!< energy of the system
		DataSet<double> corr_;			//!< correlation between neighbours
		DataSet<double> lr_corr_;		//!< long range correlation 

		Vector<unsigned int> set_ref(Parseur& P);
};
#endif
