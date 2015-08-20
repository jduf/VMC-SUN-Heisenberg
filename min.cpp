/*! @file min.cpp */

#include "VMCPSO.hpp"
#include "VMCInterpolation.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	
	VMCMinimization m(P);
	if(m.ready()){
		unsigned int i;
		switch(P.find("what",i,true)?P.get<unsigned int>(i):666){
			case 0:
				{
					VMCPSO m1(P,m);
					VMCInterpolation m2(m);
					unsigned int loop(P.get<unsigned int>("loop"));
					if(!P.locked()){
						for(unsigned int i(0);i<loop;i++){
							m1.init(true, i<3?false:true);
							m1.run();
							m1.save();

							m2.init();
							m2.run(true);
							m2.save();
						}
						//m.refine();
						m.complete_analysis(1e-5);
						m.save();
						m.plot();
						m2.plot();
					} else { std::cerr<<"min : some argument are not corretly set"<<std::endl; }
				} break;
			case 1:
				{
					unsigned int i(0);
					unsigned int j(0);
					if(P.find("E",i,false),P.find("dE",j,false)){ m.refine(P.get<double>(i),P.get<double>(j)); }
					else { m.refine(); }
					m.save();
					m.plot();
				}break;
			case 2:
				{ 
					m.set_phase_space(P); 
					m.save();
				} break;
			case 3:
				{ 
					m.plot();
				} break;
			case 4:
				{ 
					VMCPSO m1(P,m);
					VMCInterpolation m2(m);
				} break;
			default:
				{
					std::cerr<<"min : unknown option 'what', options are :"<<std::endl;
					std::cerr<<"    - complete run                       : 0"<<std::endl;
					std::cerr<<"    - refine + save + plot               : 1"<<std::endl;
					std::cerr<<"    - redefine phase space + save        : 2"<<std::endl;
					std::cerr<<"    - plot                               : 3"<<std::endl;
					std::cerr<<"    - load VMCPSO and VMCInterpolation   : 4"<<std::endl;
				}
		}
	}
}
