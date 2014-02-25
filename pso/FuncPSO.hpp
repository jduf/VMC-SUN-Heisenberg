#ifndef DEF_FUNCPSO
#define DEF_FUNCPSO

#include"PSO.hpp"
#include"Function.hpp"

class FuncPSO : public PSO{
	public:
		FuncPSO(unsigned int Nparticle, unsigned int Nparam, double cg, double cp, unsigned int maxiter);
		double f(Vector<double> x);

		void result();

	private:
		double run(Vector<double> const& x);
};

FuncPSO::FuncPSO(unsigned int Nparticle, unsigned int Nparam, double cg, double cp, unsigned int maxiter):
	PSO(Nparticle,Nparam,cg,cp,maxiter)
{}

double FuncPSO::run(Vector<double> const& x){
	Function F;
	return F.f(x);
}

void FuncPSO::result(){
	//std::cout<<"expected result ";
	//for(unsigned int i(0);i<9;i++){
		//std::cout<<Function::fsol[i]<<" ";
	//}
	//std::cout<<std::endl;
	std::cout<<"best result     "<<b_<<std::endl;
}
#endif
