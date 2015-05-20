#ifndef DEF_PSPLINE
#define DEF_PSPLINE

#include "Vector.hpp"
#include "Rand.hpp"
#include "Gnuplot.hpp"
#include "Lapack.hpp"

class PSpline{
	public:
		PSpline(unsigned int const& k);
		PSpline(unsigned int const& k, Matrix<double> const& c, Vector<double> const& y);
		void compute_weights();

		double extrapolate(Vector<double> const& x);

	private:
		std::vector<Vector<double> > c_;
		std::vector<double> y_;
		Vector<double> weights_;
		unsigned int dim_;
		unsigned int N_;
		double (PSpline::*phi_)(double const&);

		double phi1(double const& r){ return r; }
		double phi2(double const& r){ return r*r*log(r); }
		double phi3(double const& r){ return r*r*r; }
		double phi4(double const& r){ return r*r*r*r*log(r); }
		double phi5(double const& r){ return r*r*r*r*r; }
		double phi6(double const& r){ return r*r*r*r*r*r*log(r); }
};
#endif

