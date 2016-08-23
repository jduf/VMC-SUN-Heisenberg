#ifndef DEF_FIT
#define DEF_FIT

#include "Vector.hpp"
#include "cminpack.h"
#include <functional>//for std::function

class Fit{
	public:
		Fit(Vector<double> const& x, Vector<double> const& y, Vector<double> const& p, std::function<double (double, const double*)> f);
		~Fit();
		/*{Forbidden*/
		Fit() = delete;
		Fit(Fit const&) = delete;
		Fit(Fit&&) = delete;
		Fit operator=(Fit) = delete;
		/*}*/

	private:
		const double* x_;//!< pointer on a constant x
		const double* y_;//!< pointer on a constant y
		std::function<double (double, const double*)> f_;
		int m_;			//!< number of measures
		int n_;			//!< number of parameters
		int lwa_;		//!< bigger than m_*n_+5*n_+m_
		int* iwa_;		//!< work array of length n_
		double* wa_;	//!< work array of length lwa_
		double* fvec_;	//!< y_ - f_(x_,p)

		static int eval(void *data, int m, int n, const double *p, double *fvec, int iflag);
};
#endif
