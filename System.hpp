#ifndef DEF_SYSTEM
#define DEF_SYSTEM

#include "Sampling.hpp"

/*!{Class that contains the information on the state*/
/*
 * status_ = 4 : System is initialized (in System)
 * status_ = 3 : System is allowed (in System)
 * status_ = 2 : Degenerated wavefunction
 * status_ = 1 : System*D is allowed (in System*D)
 * status_ = 0 : Found an initial state (in SystemFermionic)
 */
/*}*/
class System{
	public:
		/*!Constructor*/
		System(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc);
		/*!Constructor that reads from file*/
		System(IOFiles& r);
		/*!Default destructor*/
		virtual ~System() = default;
		/*{Forbidden*/
		System(System&&) = delete;
		System& operator=(System) = delete;
		/*}*/

		/*!Returns energy*/
		Data<double> const& get_energy() const {return E_;}
		/*!Returns a ref to E_ (needed for merging)*/
		Data<double>& get_energy() {return E_;}
		/*!Returns correlation*/
		DataSet<double> const& get_corr() const {return corr_;}
		/*!Returns long range correlation*/
		DataSet<double> const& get_lr_corr() const {return lr_corr_;}
		/*!Returns the status of the system*/
		unsigned int const& get_status(){return status_;}

		/*!Sets the observables to default (0) values and initilizes binning*/
		void set_binning();
		/*!Sets the binning for E_(>=0), corr_(>=1), lr_corr(>=2)*/
		void set_observable(unsigned int const& what);
		/*!Deletes the binning for all observables*/
		void delete_binning();

		virtual void write(IOFiles& w) const;

	protected:
		/*!Copy constructor*/
		System(System const&) = default;
		/*!Default constructor*/
		System():ref_(0),N_(0),m_(0),n_(0),bc_(0){std::cout<<"System::System() : should never be called"<<std::endl;}

		Vector<unsigned int> const ref_;//!< type of system 
		unsigned int const N_;			//!< number of colors
		unsigned int const m_;			//!< number of particles per site
		unsigned int const n_;			//!< number of sites
		Vector<unsigned int> const M_;	//!< number of particles of each color 
		int const bc_;					//!< boundary condition
		unsigned int status_;			//!< status of the simulation

		Data<double> E_; 				//!< energy of the system
		DataSet<double> corr_;			//!< correlation between neighbours
		DataSet<double> lr_corr_;		//!< long range correlation 

		Matrix<unsigned int> links_;	//!< bonds <i,j> exchanged by H
};
#endif
