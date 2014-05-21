#ifndef DEF_SYSTEM
#define DEF_SYSTEM

#include "SamplingSet.hpp"
#include "Matrix.hpp"

/*!Class that contains the information on the state*/
class System{
	public:
		System(unsigned int N, unsigned int n, unsigned int m, int bc);
		System(System const& s);
		virtual ~System();

		unsigned int const& get_n() const { return n_;}
		unsigned int const& get_N() const { return N_;}
		unsigned int const& get_m() const { return m_;}
		unsigned int const& get_M() const { return M_;}
		int const& get_bc() const { return bc_;}

		System const& get_system() const { return (*this);}

		/*!Returns energy*/
		CorrelatedSamples<double> const& get_energy() const {return E_;}
		/*!Returns correlation*/
		CorrelatedSamplesSet<double> const& get_corr() const {return corr_;}
		/*!Returns long range correlation*/
		CorrelatedSamplesSet<double> const& get_long_range_corr() const {return long_range_corr_;}

		void save(IOFiles& w) const;

	protected:
		unsigned int n_;//!< number of sites
		unsigned int N_;//!< number of colors
		unsigned int m_;//!< number of particles per site
		unsigned int M_;//!< number of particles of each color 
		int bc_;		//!< boundary condition

		Matrix<unsigned int> links_;	//!< list of links

		CorrelatedSamples<double> E_;
		CorrelatedSamplesSet<double> corr_;	
		CorrelatedSamplesSet<double> long_range_corr_;

	private:
		System& operator=(System const& s);
		System();
};
#endif
