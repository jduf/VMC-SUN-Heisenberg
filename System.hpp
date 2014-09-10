#ifndef DEF_SYSTEM
#define DEF_SYSTEM

#include "Sampling.hpp"

/*!Class that contains the information on the state*/
class System{
	public:
		System(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc);
		virtual ~System(){}

		/*!Returns energy*/
		Data<double> const& get_energy() const {return E_;}
		Data<double>& get_energy() {return E_;}
		/*!Returns correlation*/
		DataSet<double> const& get_corr() const {return corr_;}
		/*!Returns long range correlation*/
		DataSet<double> const& get_lr_corr() const {return lr_corr_;}
		/*!Returns the status of the system*/
		unsigned int const& get_status(){return status_;}

		void set();

	protected:
		System(System const& s);
		System():ref_(0),N_(0),m_(0),n_(0),bc_(0){std::cout<<"System::System() : should not be called"<<std::endl;};

		/*variables that will be saved*/
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

	private:
		/*!Forbids assignment*/
		System& operator=(System const& s);
};
#endif
