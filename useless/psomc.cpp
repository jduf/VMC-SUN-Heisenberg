/*! @file psomc.cpp */

#include "VMCPSO.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);

	
	VMCMinimization m(P);
	VMCPSO s(P,m);
	std::cout<<"new independant run"<<std::endl;
	s.init(true, false);
	s.run();
	s.complete_analysis(1e-5);
	s.refine(50,0.001,10);
	s.save();
	s.plot();
	s.save_best(10);
}
