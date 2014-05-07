/*!  @file check.cpp */

#include "CreateSystem.hpp"
#include "MonteCarlo.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	CreateSystem cs(P);
	unsigned int i(0);
	if(P.search("lrc",i)){
		double param(P.get<double>("param"));
		unsigned int lrc(P.get<unsigned int>(i));
		unsigned int tmax(30);
		CreateSystem CS(cs,param);
		if(CS.use_complex()){
			MonteCarlo<std::complex<double> > sim(&CS,tmax,lrc);
			sim.check();
		} else {
			MonteCarlo<double> sim(&CS,tmax,lrc);
			sim.check();
		}
	} else {
		cs.check();
	}
}
