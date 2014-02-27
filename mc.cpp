/*!  @file mc.cpp */

#include "ParallelMonteCarlo.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int nruns(P.get<unsigned int>("nruns"));
	unsigned int tmax(P.get<unsigned int>("tmax"));
	double param(P.get<double>("param"));
	CreateSystem CS(P);
	Write w(CS.get_filename()+"-mc.jdbin");
	if( CS.use_complex() ){
		ParallelMonteCarlo<std::complex<double> > sim(CS,param,nruns,tmax);
		sim.run();
		sim.save(w);
	} else {
		ParallelMonteCarlo<double> sim(CS,param,nruns,tmax);
		sim.run();
		sim.save(w);
	}
}
