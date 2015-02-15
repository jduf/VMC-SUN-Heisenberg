#ifndef DEF_FUNCPSO
#define DEF_FUNCPSO

#include"PSO.hpp"

class FuncPSO : public PSO{
	public:
		FuncPSO(unsigned int Nparticle, unsigned int Nparam, double cg, double cp, unsigned int maxiter);
		double f(Vector<double> x);

		void result();

	private:
		double f(Vector<double> const& x);
};

FuncPSO::FuncPSO(unsigned int Nparticle, unsigned int Nparam, double cg, double cp, unsigned int maxiter):
	PSO(Nparticle,Nparam,cg,cp,maxiter)
{}

double FuncPSO::f(Vector<double> const& x){
	return (1-x(0)*x(0))*(1-x(0)*x(0))+(3-x(1)*x(1))*(3-x(1)*x(1));
}

void FuncPSO::result(){}
#endif
