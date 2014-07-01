#ifndef DEF_SYSTEM
#define DEF_SYSTEM

#include "Sampling.hpp"

/*!Class that contains the information on the state*/
class System{
	public:
		System(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc);
		virtual ~System(){
	std::cout<<"destroy System"<<std::endl;
		}

		unsigned int const& get_N() const { return N_;}
		unsigned int const& get_m() const { return m_;}
		unsigned int const& get_n() const { return n_;}
		int const& get_bc() const { return bc_;}

		System const* get_system() const { return this;}

		/*!Returns energy*/
		Data<double> const& get_energy() const {return E_;}
		/*!Returns correlation*/
		DataSet<double> const& get_corr() const {return corr_;}
		/*!Returns long range correlation*/
		DataSet<double> const& get_long_range_corr() const {return long_range_corr_;}

		void save(IOFiles& w) const;

	protected:
		System(System const& s);
		System():ref_(0),N_(0),m_(0),n_(0),bc_(0){std::cout<<"system default"<<std::endl;};

		/*variables that will be saved*/
		Vector<unsigned int> const ref_;
		unsigned int const N_;//!< number of colors
		unsigned int const m_;//!< number of particles per site
		unsigned int const n_;//!< number of sites
		Vector<unsigned int> M_;//!< number of particles of each color 
		int const bc_;//!< boundary condition

		Data<double> E_;
		DataSet<double> corr_;
		DataSet<double> long_range_corr_;

		Matrix<unsigned int> links_;//!< list of links

	private:
		System& operator=(System const& s);
};
#endif
