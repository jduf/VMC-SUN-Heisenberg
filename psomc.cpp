/*!  @file psomc.cpp */

#include "PSOMonteCarlo.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	
	unsigned int Nfreedom(P.get<unsigned int>("Nfreedom"));

	PSOMonteCarlo s(&P);
	unsigned int i(0);
	if(P.find("load",i,false)){
		IOFiles r(P.get<std::string>(i),false);
		s.read(r); 
	}

	for(unsigned int i(0);i<Nfreedom;i++){
		Optimization::set_limit(i,-1.2,-0.8);
	}
	s.init(100);
	s.run();
	s.complete_analysis(1e-5);
	s.print();
	s.plot();
	if(P.find("save",i,false)){
		IOFiles w(P.get<std::string>(i),true);
		s.write(w);
	}
}
