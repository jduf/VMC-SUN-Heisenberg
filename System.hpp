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
		/*!Returns correlation*/
		DataSet<double> const& get_corr() const {return corr_;}
		/*!Returns long range correlation*/
		DataSet<double> const& get_long_range_corr() const {return long_range_corr_;}
		/*!Returns the status of the system*/
		unsigned int const& get_status(){return status_;}

		void save(IOFiles& w) const;

		void set();

	protected:
		System(System const& s);
		System():ref_(0),N_(0),m_(0),n_(0),bc_(0){std::cout<<"system default"<<std::endl;};

		/*variables that will be saved*/
		Vector<unsigned int> const ref_;
		unsigned int const N_;	//!< number of colors
		unsigned int const m_;	//!< number of particles per site
		unsigned int const n_;	//!< number of sites
		Vector<unsigned int> const M_;//!< number of particles of each color 
		int const bc_;			//!< boundary condition
		unsigned int status_;	//!< 

		Data<double> E_;
		DataSet<double> corr_;
		DataSet<double> long_range_corr_;

		Matrix<unsigned int> links_;//!< list of links

	private:
		/*!Forbids assignment*/
		System& operator=(System const& s);
};
#endif
