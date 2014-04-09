#ifndef DEF_FIT
#define DEF_FIT

#include "Vector.hpp"
#include "levmar.h"

class Fit{
	public:
		Fit(Vector<double> const& x, Vector<double> const& y, double (*f)(double,double*), Vector<double>& p);

		Vector<double> fx() const;

	private:
		static void func(double *p, double *y, int m, int n, void *adata);
		double operator()(unsigned int i, double* p) const { return f_(x_(i),p); }

		double (*f_)(double,double*);
		Vector<double> x_;
		Vector<double> p_;
		int ret_;
};
#endif
