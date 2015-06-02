/*! @file psomc.cpp */

#include "VMCPSO.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);

	
	Minimization m(P);
	VMCPSO s(P,m);
	//for(unsigned int i(0);i<Optimization::get_Nfreedom();i++){
		//Optimization::set_limit(i,-2.0,2.0);
	//}
	//Optimization::set_limit(0,-1.0,-0.5);
	//Optimization::set_limit(1,0.7,1.3);
	//Optimization::set_limit(2,-1.3,-0.7);
	s.set_ps(0,P.get<std::vector<double> >("t1"));
	s.set_ps(1,P.get<std::vector<double> >("t2"));
	s.set_ps(2,P.get<std::vector<double> >("t3"));
	std::cout<<"new independant run"<<std::endl;
	s.init(true, false);
	s.run();
	s.complete_analysis(1e-5);
	s.refine(50,0.001,10);
	s.save();
	s.print();
	s.plot();
	s.save_best(10);
}
