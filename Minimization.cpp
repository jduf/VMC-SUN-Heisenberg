#include "Minimization.hpp"

Minimization::Minimization(Parseur& P):
	nthreads(1),
	t_max(P.get<unsigned int>("t_max")),
	CS_(P)
{
	P.get("nthreads",nthreads);
}

Minimization::~Minimization(){}

void Minimization::min(double xmax){
	double x0(0.0);
	double f0(fx(x0));

	double x1(xmax);
	double f1(fx(x1));

	double x((x0+x1)/2.0);
	double f(fx(x));

	double xt0(0.0);
	double ft0(0.0);

	double xt1(0.0);
	double ft1(0.0);
	for(unsigned int i(0);i<5;i++){
		xt0 = (x0+x)/2.0;
		xt1 = (x1+x)/2.0;
		ft0 = fx(xt0);
		ft1 = fx(xt1);
		if((f0-ft0)*(ft0-f)>0 && (f-ft1)*(ft1-f1)>0 ){
			f0=ft0;
			x0=xt0;
			f1=ft1;
			x1=xt1;
		} 
		if((f0-ft0)*(ft0-f)<0 && (f-ft1)*(ft1-f1)>0 ){
			f1=f;
			x1=x;
			f=ft0;
			x=xt0;
		} 
		if((f0-ft0)*(ft0-f)>0 && (f-ft1)*(ft1-f1)<0 ){
			f0=f;
			x0=x;
			f=ft1;
			x=xt1;
		} 
	}
}

double Minimization::fx(double delta){
	double energy(0.0);
#pragma omp parallel num_threads(nthreads)
	{
		CreateSystem cs(CS_,delta);
		MonteCarlo<double> sim(cs,t_max);
		sim.run(2);
		std::cout<<delta<<" "<<sim.get_energy()<<" "<<sim.get_error()<<std::endl;
		energy+=sim.get_energy();
	}
	std::cout<<delta<<" "<<energy / nthreads<<std::endl;
	return energy / nthreads; 
}
