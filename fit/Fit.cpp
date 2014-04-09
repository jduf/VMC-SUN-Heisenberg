#include "Fit.hpp"

Fit::Fit(Vector<double> const& x, Vector<double> const& y, double (*f)(double, double*), Vector<double>& p):
	f_(f),
	x_(x),
	p_(p),
	ret_(0)
{
	double info[10];
	ret_ = dlevmar_dif(func,p_.ptr(),y.ptr(),p_.size(),y.size(),1e2,NULL,info,NULL,NULL,this);
	p = p_;
}

void Fit::func(double *p, double *y, int m, int n, void *adata){
	Fit* self = static_cast<Fit*>(adata);
	for(int i(0);i<n;i++){ y[i] = (*self)(i,p); }
}

Vector<double> Fit::fx() const {
	Vector<double> y(x_.size());
	for(unsigned int i(0);i<y.size();i++){
		y(i) = f_(x_(i),p_.ptr());
	}
	return y;
}
