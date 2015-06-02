/*! @file min.cpp */

#include "VMCPSO.hpp"
#include "VMCSpline.hpp"


int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	
	VMCPSO m1(P);
	VMCSpline m2(P);
	m1.set_ps(0,P.get<std::vector<double> >("t1"));
	m1.set_ps(1,P.get<std::vector<double> >("t2"));
	m1.set_ps(2,P.get<std::vector<double> >("t3"));
	m2.set_ps(0,P.get<std::vector<double> >("t1"));
	m2.set_ps(1,P.get<std::vector<double> >("t2"));
	m2.set_ps(2,P.get<std::vector<double> >("t3"));

	for(unsigned int i(0);i<3;i++){
		if(i){m1.move(&m2); }
		m1.init(true, false);
		m1.run();
		//m1.save();
		//m1.print();
		//m1.plot();

		m2.move(&m1);
		m2.init();
		m2.run(1);
		m2.plot();
	}

	m2.complete_analysis(1e-5);
	m2.refine(50,0.001,10);
	m2.print();
	m2.save();
}

//int main(int argc, char* argv[]){
	//Parseur P(argc,argv);
	//
	//VMCSpline m2(P);
	//m2.set_ps(0,P.get<std::vector<double> >("t1"));
	//m2.set_ps(1,P.get<std::vector<double> >("t2"));
	//m2.set_ps(2,P.get<std::vector<double> >("t3"));
//
	//m2.init();
	//m2.run(1);
	//m2.plot();
	//m2.print();
//}
