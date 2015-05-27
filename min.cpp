/*! @file min.cpp */

#include "VMCPSO.hpp"
#include "VMCSpline.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	
	//VMCPSO pso(P);
	//for(unsigned int i(0);i<Optimization::get_Nfreedom();i++){
		//Optimization::set_limit(i,-2.0,2.0);
	//}
	////Optimization::set_limit(0,-1.0,-0.5);
	////Optimization::set_limit(1,0.7,1.3);
	////Optimization::set_limit(2,-1.3,-0.7);
	//std::cout<<"new independant run"<<std::endl;
	//pso.init(true, false);
	//pso.run(0.8);
	//pso.complete_analysis(1e-5);
	////pso.refine(50,0.001,10);
	//pso.save();
	//pso.print();
	//pso.plot();

	VMCSpline spline(P);
	spline.set_x(0,P.get<std::vector<double> >("t1"));
	spline.set_x(1,P.get<std::vector<double> >("t2"));

	//spline.move(&pso);
	spline.run();
	spline.plot();
}
