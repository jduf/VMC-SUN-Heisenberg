/*!  @file mc.cpp */

#include "ParallelMonteCarlo.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int nruns(P.get<unsigned int>("nruns"));
	unsigned int tmax(P.get<unsigned int>("tmax"));
	Vector<double> param(P.get<Vector<double> >("param"));
	CreateSystem cs(P);
	for(unsigned int i(0);i<param.size();i++){
		CreateSystem CS(cs,param(i));
		Write w(CS.get_filename()+"-mc.jdbin");
		if( CS.use_complex() ){
			ParallelMonteCarlo<std::complex<double> > sim(&CS,nruns,tmax);
			sim.run();
			sim.save(w);
		} else {
			ParallelMonteCarlo<double> sim(&CS,nruns,tmax);
			sim.run();
			sim.save(w);
		}
	}
}
