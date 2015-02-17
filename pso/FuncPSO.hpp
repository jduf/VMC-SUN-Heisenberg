#ifndef DEF_FUNCPSO
#define DEF_FUNCPSO

#include "PSO.hpp"
#include "List.hpp"
#include "Vector.hpp"

class Measure{
	public:
		Measure(Vector<double> const& x, double const& fx):
			x_(x), fx_(fx), N_(0){};

		void print(std::ostream& flux) const
		{ flux<<x_<<" : "<<fx_<<std::endl; }

		Vector<double> x_;
		double fx_;
		unsigned int N_;
};

bool func(Measure* a, Measure* b) { 
	unsigned int i(0);
	while(i<a->x_.size()){
		if(a->x_(i) > b->x_(i)){ return false; }
		if(a->x_(i) < b->x_(i)){ return true; }
		if(a->x_(i)== b->x_(i)){ i++; }
	}
	return false;
}

std::ostream& operator<<(std::ostream& flux, Measure const& m){
	m.print(flux);
	return flux;
}

class FuncPSO : public PSO {
	public:
		FuncPSO(unsigned int Nparticle, unsigned int Nparam, double cg, double cp, unsigned int maxiter);
		double f(Vector<double> x);

		void result();

		double f(Vector<double> const& x);
		List<Measure> m_;
};

FuncPSO::FuncPSO(unsigned int Nparticle, unsigned int Nparam, double cg, double cp, unsigned int maxiter):
	PSO(Nparticle,Nparam,cg,cp,maxiter)
{}

double FuncPSO::f(Vector<double> const& x){
	double fx((1-x(0)*x(0))*(1-x(0)*x(0))+(3-x(1)*x(1))*(3-x(1)*x(1)));
#pragma omp critical
	{
		m_.add_sort(new Measure(x,fx),func);
	}
	return fx;
}

void FuncPSO::result(){
	PSO_print();
	std::cout<<m_<<std::endl;
}
#endif
