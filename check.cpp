/*!  @file check.cpp */

#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	CreateSystem cs(&P);
	unsigned int i(0);
	switch(P.find("what",i)?P.get<unsigned int>(i):0){
		case 1:/*call CreateSystem::init and check*/
			{ 
				cs.init();
				if(cs.get_status()==2){
					cs.check();
				}
			} break;
		case 2:/*call CreateSystem::init create and check*/
			{ 
				cs.init();
				if(cs.get_status()==2){
					cs.create();
					cs.check();
				}
			} break;
		case 3:/*call CreateSystem::init create and run MonteCarlo*/
			{
				unsigned int tmax(P.find("tmax",i,false)?P.get<unsigned int>(i):10);
				cs.init();
				if(cs.get_status()==2){
					cs.create();
					MCSystem* S(NULL);
					if(cs.use_complex()){
						S = new SystemFermionic<std::complex<double> >(*dynamic_cast<const Fermionic<std::complex<double> >*>(cs.get_system())); 
					} else {
						S = new SystemFermionic<double>(*dynamic_cast<const Fermionic<double>*>(cs.get_system())); 
					}
					MonteCarlo sim(S,tmax);
					sim.run();
					std::cout<<S->get_energy()<<std::endl;
					delete S;
				}
			} break;
		default:
			{
				std::cerr<<"check : unknown option 'what', options are :"<<std::endl;
				std::cerr<<"      - init + check          : 1"<<std::endl;
				std::cerr<<"      - init + create + check : 2"<<std::endl;
				std::cerr<<"      - init + create + run   : 3"<<std::endl;
			}
	} 
}
