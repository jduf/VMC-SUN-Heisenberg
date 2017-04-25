#include "Minimization.hpp"

Minimization::Minimization(Parseur& P):
	CS_(P),
	locked_(false),
	tmax_(P.get<unsigned int>("tmax")),
	nruns_(P.get<unsigned int>("nruns"))
{
	std::cout<<"need to define the correct path"<<std::endl;
	if(P.status()){
		std::cerr<<"Minimization::Minimization(Parseur& P) : An argument is missing, the minimization won't be lunched"<<std::endl;
		locked_=true;
	}
}

Minimization::~Minimization(){}

void Minimization::min(double xmax){
	if(!locked_){
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
		unsigned int i(0);
		unsigned int cond;
		do{
			cond = 0;
			xt0 = (x0+x)/2.0;
			xt1 = (x1+x)/2.0;
			ft0 = fx(xt0);
			ft1 = fx(xt1);
			if((f0-ft0)*(ft0-f)>0 && (f-ft1)*(ft1-f1)>0 ){
				f0=ft0;
				x0=xt0;
				f1=ft1;
				x1=xt1;
				cond++;
			}
			if((f0-ft0)*(ft0-f)<0 && (f-ft1)*(ft1-f1)>0 ){
				f1=f;
				x1=x;
				f=ft0;
				x=xt0;
				cond++;
			}
			if((f0-ft0)*(ft0-f)>0 && (f-ft1)*(ft1-f1)<0 ){
				f0=f;
				x0=x;
				f=ft1;
				x=xt1;
				cond++;
			}
		}while (++i<10 && cond==1);
	}
}

double Minimization::fx(double param){
	results_file_("param",param);
	CreateSystem CS(CS_,param);
	if( CS.use_complex() ){
		ParallelMonteCarlo<std::complex<double> > sim(&CS,".",nruns_,tmax_);
		sim.run(1);
		return sim.get_energy();
	} else {
		ParallelMonteCarlo<double> sim(&CS,".",nruns_,tmax_);
		sim.run(1);
		return sim.get_energy();
	}
}
