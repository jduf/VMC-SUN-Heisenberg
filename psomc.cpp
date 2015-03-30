/*!  @file psomc.cpp */

#include "PSOFermionic.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	
	unsigned int Nfreedom(P.get<unsigned int>("Nfreedom"));

	PSOFermionic s(&P);

	for(unsigned int i(0);i<Nfreedom;i++){
		Optimization::set_limit(i,0.9,1.1);
	}
	s.init(100);
	s.run();
	s.print();
	s.plot();
}
