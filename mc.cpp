/*!  @file mc.cpp */

#include "MonteCarlo.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int nthreads(1);
	P.get("nthreads",nthreads);
	CreateSystem CS(P);
	double param(P.get<double>("param"));

	if( !CS.use_complex() ){
#pragma omp parallel num_threads(nthreads)
		{
			CreateSystem cs(CS,param);
			MonteCarlo<double> sim(cs);
			sim.run(1,1e6);
			std::cout<<param<<" "<<sim.get_energy()<<" "<<sim.get_error()<<std::endl;
			cs.save(sim.get_energy(),sim.get_error(),sim.get_correlation());
		}
	} else {
#pragma omp parallel num_threads(nthreads)
		{
			CreateSystem cs(CS,param);
			MonteCarlo<std::complex<double> > sim(cs);
			sim.run(1,1e6);
			std::cout<<param<<" "<<sim.get_energy()<<" "<<sim.get_error()<<std::endl;
			cs.save(sim.get_energy(),sim.get_error(),sim.get_correlation());
		}
	}
}
