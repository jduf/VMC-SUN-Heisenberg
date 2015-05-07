/*!  @file psomc.cpp */

#include "PSOMonteCarlo.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	
	unsigned int Nfreedom(P.get<unsigned int>("Nfreedom"));
	unsigned int i(0);
	IOFiles* in(P.find("load",i,false)?(new IOFiles(P.get<std::string>(i),false)):NULL);
	PSOMonteCarlo s(P,in);
	if(in){ 
		delete in;
		in = NULL;
	}
	for(unsigned int i(0);i<Nfreedom;i++){
		Optimization::set_limit(i,-2.0,2.0);
	}
	Optimization::set_limit(0,-1.0,-0.5);
	Optimization::set_limit(1,0.7,1.3);
	Optimization::set_limit(2,-1.3,-0.7);

	s.create_particle_history(true); 
	s.init(100);
	s.run();
	s.complete_analysis(1e-5);
	s.refine(500,0.01,10);
	s.save();
	s.print();
	s.plot();
}
