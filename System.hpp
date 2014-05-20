#ifndef DEF_SYSTEM
#define DEF_SYSTEM

#include "SamplingSet.hpp"
#include "Matrix.hpp"

/*!Class that contains the information on the state*/
class System{
	public:
		System(unsigned int N, unsigned int n, unsigned int m, int bc);
		System(System const& S);
		virtual ~System();

		/*!Returns energy*/
		CorrelatedSamples<double> const& get_energy() const {return E_;}
		/*!Returns correlation*/
		CorrelatedSamplesSet<double> const& get_corr() const {return corr_;}
		/*!Returns long range correlation*/
		CorrelatedSamplesSet<double> const& get_long_range_corr() const {return long_range_corr_;}

		void save(IOFiles& w) const;

	protected:
		unsigned int const n_;		//!< number of sites
		unsigned int const N_;		//!< number of colors
		unsigned int const m_;		//!< number of particles per site
		unsigned int const M_;		//!< number of particles of each color 
		int bc_;					//!< boundary condition

		Matrix<unsigned int> links_;	//!< list of links

		CorrelatedSamples<double> E_;
		CorrelatedSamplesSet<double> corr_;	
		CorrelatedSamplesSet<double> long_range_corr_;
};
#endif
