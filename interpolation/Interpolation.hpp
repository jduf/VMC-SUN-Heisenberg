#ifndef DEF_INTERPOLATION
#define DEF_INTERPOLATION

#include "Lapack.hpp"
#include "omp.h"

/*!If two data <c,y> are identical the matrix to invert won't be invertible*/
class Interpolation{
	public:
		Interpolation(unsigned int const& dim);

		bool select_basis_function(unsigned int const& basis);
		void add_data(Vector<double> const& ci, double const& yi);
		void add_data(unsigned int const& i, Vector<double> const& c, double const& y);

		void set_data();
		void compute_weights();
		bool compute_weights(double const& dx, unsigned int const& n);

		double extrapolate(Vector<double> const& x) const;

	private:
		std::vector<Vector<double> > c_;
		std::vector<double> y_;
		Vector<double> weights_;
		unsigned int basis_;
		unsigned int dim_;
		unsigned int N_;
		double support_;
		double (Interpolation::*phi_)(double const&) const;

		double phi1(double const& r) const { return r; }
		double phi2(double const& r) const { return r*(r<1?log(pow(r,r)):r*log(r)); }
		double phi3(double const& r) const { return r*r*r; }
		double phi4(double const& r) const { return r*r*r*(r<1?log(pow(r,r)):r*log(r)); }
		double phi5(double const& r) const { return r*r*r*r*r; }
		double phi6(double const& r) const { return r*r*r*r*r*(r<1?log(pow(r,r)):r*log(r)); }
		double phi7(double const& r) const { return (r<1.0?(1-r)*(1-r):0.0); }
		double phi8(double const& r) const { return (r<1.0?(1-r)*(1-r)*(1-r)*(1-r)*(4*r+1):0.0); }
		double phi9(double const& r) const { return (r<1.0?(1-r)*(1-r)*(1-r)*(1-r)*(1-r)*(1-r)*(32*r*r*r+25*r*r+8*r+1):0.0); }
};
#endif
