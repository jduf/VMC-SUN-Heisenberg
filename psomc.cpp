/*!  @file psomc.cpp */

#include "PSOFermionic.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	
	unsigned int Nfreedom(1);

	PSOFermionic s(P);

	for(unsigned int i(0);i<Nfreedom;i++){
		s.PSO_set_limit(i,-0.3,0.3);
	}
	s.PSO_init();
	s.PSO_run();
	s.PSO_print();
}
