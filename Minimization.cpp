#include "Minimization.hpp"

Minimization::Minimization(Parseur& P):
	CS_(P),
	tmax_(P.get<unsigned int>("tmax")),
	nruns_(P.get<unsigned int>("nruns")),
	results_file_(CS_.get_filename()+"-min.jdbin")
{
	results_file_("Created by the Minimization class",1);
	results_file_("nruns (number of simulations runned)",nruns_);
	RST rst_param;
	rst_param.title("Input","-");
	results_file_.add_to_header(rst_param.get());
	CS_.save(results_file_);
	RST rst_results;
	rst_results.title("Results","-");
	results_file_.add_to_header(rst_results.get());
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

double Minimization::fx(double param){
	results_file_("param",param);
	if( CS_.use_complex() ){
		ParallelMonteCarlo<std::complex<double> > sim(CS_,param,nruns_,tmax_);
		sim.run();
		sim.save(results_file_);
		return sim.get_energy();
	} else {
		ParallelMonteCarlo<double> sim(CS_,param,nruns_,tmax_);
		sim.run();
		sim.save(results_file_);
		return sim.get_energy();
	}
}
