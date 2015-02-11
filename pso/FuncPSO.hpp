#ifndef DEF_FUNCPSO
#define DEF_FUNCPSO

#include"PSO.hpp"

class FuncPSO : public PSO{
	public:
		FuncPSO(unsigned int Nparticle, unsigned int Nparam, double cg, double cp, unsigned int maxiter);
		double f(Vector<double> x);

		void result();

	private:
		double evaluate(Vector<double> const& x);
		void set_to_grid(unsigned int i);
};

FuncPSO::FuncPSO(unsigned int Nparticle, unsigned int Nparam, double cg, double cp, unsigned int maxiter):
	PSO(Nparticle,Nparam,cg,cp,maxiter)
{
	for(unsigned int i(0);i<Nfreedom_;i++){ PSO_set_limit(i,-2,2); }
}

double FuncPSO::evaluate(Vector<double> const& x){
	return (1-x(0)*x(0))*(1-x(0)*x(0))+(3-x(1)*x(1))*(3-x(1)*x(1));
}

void FuncPSO::result(){
	std::cout<<"best result     "<<pfb_[bbee_]<<" : "<<pb_[bbee_]<<std::endl;
}

void FuncPSO::set_to_grid(unsigned int i){
	unsigned int n;
	double dx(0.01);
	for(unsigned int j(0);j<Nfreedom_;j++){
		n=0;
		if(std::abs(px_[i](j))<dx/2){ n=1; }
		if(std::abs(px_[i](j)-min_(j))<dx/2){ n=2; }
		if(std::abs(px_[i](j)-max_(j))<dx/2){ n=3; }
		switch(n){
			case 0:{ px_[i](j) = std::round(std::abs(px_[i](j)/dx))*dx; }break;
			case 1:{ px_[i](j) = 0; }break;
			case 2:{ px_[i](j) = min_(j); }break;
			case 3:{ px_[i](j) = max_(j); }break;
		}
	}
}
#endif
