#ifndef DEF_SYSTEM
#define DEF_SYSTEM

#include "Sampling.hpp"
#include "Matrix.hpp"

/*!Class that contains the information on the state*/
class System{
	public:
		System(unsigned int const& N, unsigned int const& n, unsigned int const& m, int const& bc, Vector<unsigned int> const& ref);
		virtual ~System(){}

		unsigned int const& get_n() const { return n_;}
		unsigned int const& get_N() const { return N_;}
		unsigned int const& get_m() const { return m_;}
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
		System():n_(0),N_(0),m_(0),ref_(0){std::cout<<"system default"<<std::endl;};

		/*variables that will be saved*/
		unsigned int const n_;//!< number of sites
		unsigned int const N_;//!< number of colors
		unsigned int const m_;//!< number of particles per site
		Vector<unsigned int> const ref_;
		Vector<unsigned int> M_;//!< number of particles of each color 
		int bc_;//!< boundary condition

		Data<double> E_;
		DataSet<double> corr_;
		DataSet<double> long_range_corr_;

		Matrix<unsigned int> links_;//!< list of links

	private:
		System& operator=(System const& s);
};
#endif
