/*!  @file psomc.cpp */

#include "VMCPSO.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	
	unsigned int Nfreedom(P.get<unsigned int>("Nfreedom"));
	VMCPSO s(&P);

	for(unsigned int i(0);i<Nfreedom;i++){
		Particle::set_limit(i,-3.0,3.0);
	}
	//Optimization::set_limit(0,-0.5,0.5);
	//Optimization::set_limit(1,.0,3.0);
	//Optimization::set_limit(2,.0,3.0);

	unsigned int i(0);
	if(P.find("load",i,false)){
		IOFiles r(P.get<std::string>(i),false);
		s.read(r,false); 
	}
	//s.init(100);
	//s.run();
	//s.complete_analysis(1e-5);
	s.refine(30,0.001,10);
	//s.print();
	s.plot();
	if(P.find("save",i,false)){
		IOFiles w(P.get<std::string>(i),true);
		s.write(w);
	}
}
