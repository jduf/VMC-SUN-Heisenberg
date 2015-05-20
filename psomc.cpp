/*!  @file psomc.cpp */

#include "PSOMonteCarlo.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	
	PSOMonteCarlo s(P);
	for(unsigned int i(0);i<Optimization::get_Nfreedom();i++){
		Optimization::set_limit(i,-2.0,2.0);
	}
	//Optimization::set_limit(0,-1.0,-0.5);
	//Optimization::set_limit(1,0.7,1.3);
	//Optimization::set_limit(2,-1.3,-0.7);
	double tol(0.8);
	for(unsigned int i(0);i<10;i++){
		std::cout<<"new independant run"<<std::endl;
		s.init(true, false);
		if(s.run(tol)){  tol *= 2./3.; }
		else { tol *= 3./2.; }
		s.complete_analysis(1e-5);
		s.refine(50,0.001,10);
		s.save();
		s.print();
		s.plot();
	}
	tol=1.0;
	for(unsigned int i(0);i<3;i++){
		std::cout<<"new independant run"<<std::endl;
		s.init(true, true);
		if(s.run(tol)){  tol /= 2.0; }
		else { tol *= 2.0; }
		s.complete_analysis(1e-5);
		s.refine(50,0.001,10);
		s.save();
		s.print();
		s.plot();
	}
	s.save(10);
}
