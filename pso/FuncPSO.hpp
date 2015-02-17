#ifndef DEF_FUNCPSO
#define DEF_FUNCPSO

#include "PSO.hpp"
#include "List.hpp"
#include <unistd.h>

class Measure{
	public:
		Measure(Vector<double> const& x; double const& fx):
			x_(x),
			fx_(fx),
			N_(0)
	{};

	private:
		Vector<double> x_;
		double fx_;
		unsigned int N_;
};

class FuncPSO : public PSO{
	public:
		FuncPSO(unsigned int Nparticle, unsigned int Nparam, double cg, double cp, unsigned int maxiter);
		double f(Vector<double> x);

		void result();

	private:
		double f(Vector<double> const& x);
		Rand<double> r_;
		List<Measure> m_;
};

FuncPSO::FuncPSO(unsigned int Nparticle, unsigned int Nparam, double cg, double cp, unsigned int maxiter):
	PSO(Nparticle,Nparam,cg,cp,maxiter),
	r_(500000,600000.0)
{}

double FuncPSO::f(Vector<double> const& x){
	double fx((1-x(0)*x(0))*(1-x(0)*x(0))+(3-x(1)*x(1))*(3-x(1)*x(1)));
	m_.append(Measure(x,fx));
	//usleep(r_.get());
	return fx;
}

void FuncPSO::result(){
	PSO_print();
}
#endif
