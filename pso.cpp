/*!  @file pso.cpp */

#include "PSOFermionic.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	
	unsigned int Nfreedom(1);

	PSOFermionic s(P,Nfreedom);

	for(unsigned int i(0);i<Nfreedom;i++){
		s.PSO_set_limit(i,0.8,1.2);
	}
	s.PSO_init();
	s.PSO_run(false); /*false because each run can vary in time*/
}
