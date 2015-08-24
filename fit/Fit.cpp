#include "Fit.hpp"

Fit::Fit(Vector<double> const& x, Vector<double> const& y, Vector<double> const& p, std::function<double (double, const double*)> f):
	x_(x.ptr()),
	y_(y.ptr()),
	f_(f),
	m_(x.size()),
	n_(p.size()),
	lwa_((m_+5)*n_+m_),
	iwa_(new int[n_]),
	wa_(new double[lwa_]),
	fvec_(new double[m_])
{
	if(m_>n_){
		double tol(sqrt(dpmpar(1)));
		lmdif1(eval, this, m_, n_, p.ptr(), fvec_, tol, iwa_, wa_, lwa_);
	} else {
		std::cerr<<__PRETTY_FUNCTION__<<" : the number of measures must be bigger than the number of parameters."<<std::endl;
	}
}

Fit::~Fit(){
	delete[] iwa_;
	delete[] wa_;
	delete[] fvec_;
}

int Fit::eval(void *data, int m, int n, const double *p, double *fvec, int iflag) {
	(void)(iflag);
	(void)(n);

	std::function<double (double, const double*)> f(((Fit*)data)->f_);
	const double* x(((Fit*)data)->x_);
	const double* y(((Fit*)data)->y_);

	for(int i(0);i<m;++i){ fvec[i] = y[i]- f(x[i],p); }
	return 0;
}
