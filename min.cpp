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

					} else { std::cerr<<__PRETTY_FUNCTION__<<" : some argument are not corretly set"<<std::endl; }
				}break;
			case 1:
				{
					unsigned int i(0);
					unsigned int j(0);
					if(P.find("E",i,false),P.find("dE",j,false)){ m.refine(P.get<double>(i),P.get<double>(j)); }
					else { m.refine(); }
					m.save();
				}break;
			case 2:
				{ 
					m.set_phase_space(P); 
					m.save();
				}break;
			case 3:
				{ 
					m.find_and_run_minima(10);
					m.save();
				}break;
			case 4:
				{ 
					IOFiles out("test.jdbin",true);
					m.find_save_and_plot_minima(5,out);
					VMCInterpolation m2(m);
					m2.plot();
				} break;
			case 5:
				{
					VMCPSO m1(P,m);
					VMCInterpolation m2(m);
				}break;
			default:
				{
					std::cerr<<__PRETTY_FUNCTION__<<" : unknown option 'what', options are :"<<std::endl;
					std::cerr<<"    - complete run                       : 0"<<std::endl;
					std::cerr<<"    - refine + save                      : 1"<<std::endl;
					std::cerr<<"    - redefine phase space + save        : 2"<<std::endl;
					std::cerr<<"    - find and run minima + save         : 3"<<std::endl;
					std::cerr<<"    - plot                               : 4"<<std::endl;
					std::cerr<<"    - load VMCPSO and VMCInterpolation   : 5"<<std::endl;
				}
		}
	}
}
