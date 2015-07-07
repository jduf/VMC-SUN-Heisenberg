/*! @file min.cpp */

#include "VMCPSO.hpp"
#include "VMCSpline.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	
	VMCMinimization m(P);
	VMCPSO m1(P,m);
	VMCSpline m2(m);
	for(unsigned int i(0);i<1;i++){
		m1.init(true, false);
		m1.run();
		m1.save();

		m2.init();
		m2.run(0);
		m2.plot();
		m2.complete_analysis(1e-5);
		m2.save();
	}
	//VMCSpline m2(m);
	//m2.init();
	//m2.run(0);
	//m2.plot();
	//m2.complete_analysis(1e-5);
}
