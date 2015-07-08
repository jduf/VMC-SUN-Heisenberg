/*! @file min.cpp */

#include "VMCPSO.hpp"
#include "VMCSpline.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int i;
	unsigned int loop(P.find("loop",i,false)?P.get<unsigned int>(i):0);
	
	VMCMinimization m(P);
	if(loop){
		VMCPSO m1(P,m);
		VMCSpline m2(m);
		if(!P.status()){
			for(unsigned int i(0);i<loop;i++){
				m1.init(true, false);
				m1.run();
				m1.save();

				m2.init();
				m2.run(0);
				m2.plot();
				m2.complete_analysis(1e-5);
				m2.save();
			}
		} else { std::cerr<<"min : some argument are not corretly set"<<std::endl; }
	} else {
		if(P.find("load",i,true)){
			VMCSpline m1(m);

			m1.init();
			m1.run(0);
			m.plot();
		}
	}
}
