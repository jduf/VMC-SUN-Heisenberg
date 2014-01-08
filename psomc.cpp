/*!  @file pso.cpp */

#include "Container.hpp"
#include "Parseur.hpp"
#include "MonteCarlo.hpp"
#include "PSOMonteCarlo.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	
	PSOMonteCarlo s(P);
	for(unsigned int i(0);i<P.get<unsigned int>("Nfreedom");i++){
		s.PSO_set_limit(i,-0.5,0.9);
	}
	s.PSO_init();
	s.PSO_run(false);
	s.save();
}


