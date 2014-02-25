/*!  @file pso.cpp */

#include "PSOFermionic.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	
	unsigned int Nfreedom(2);

	PSOFermionic s(P,Nfreedom,1e6);

	for(unsigned int i(0);i<Nfreedom;i++){
		s.PSO_set_limit(i,0.9,1.1);
	}
	s.PSO_init();
	s.PSO_run(false); /*false because each run can vary in time*/
	s.PSO_print();
	s.PSO_save("test.jdbin");
	PSOFermionic s2(P,Nfreedom,1e7);
	for(unsigned int i(0);i<Nfreedom;i++){
		s2.PSO_set_limit(i,0.9,1.1);
	}
	s2.PSO_load("test.jdbin");
	s2.PSO_run(false);
	s2.PSO_print();

}
