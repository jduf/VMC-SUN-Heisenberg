/*!  @file check.cpp */

#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	std::cout<<"############# Init System #################"<<std::endl;
	System s(P);
	std::cout<<"############# Init CreateSystem ###########"<<std::endl;
	CreateSystem cs(&s);
	unsigned int i(0);
	switch(P.find("what",i)?P.get<unsigned int>(i):0){
		case 0:/*call CreateSystem::init*/
			{ 
				cs.set_param(&P,NULL);
				std::cout<<"############# Init GenericSystem ##########"<<std::endl;
				cs.construct_GenericSystem(NULL,NULL); 
			} break;
		case 1:/*call CreateSystem::init and check*/
			{ 
				cs.set_param(&P,NULL);
				std::cout<<"############# Init GenericSystem ##########"<<std::endl;
				cs.construct_GenericSystem(NULL,NULL);
				std::cout<<"############# Check #######################"<<std::endl;
				if(cs.get_status()==2){ cs.check(); }
			} break;
		case 2:/*call CreateSystem::init create and check*/
			{ 
				cs.set_param(&P,NULL);
				std::cout<<"############# Init GenericSystem ##########"<<std::endl;
				cs.construct_GenericSystem(NULL,NULL);
				if(cs.get_status()==2){
					std::cout<<"############# Create GenericSystem ########"<<std::endl;
					cs.create();
					std::cout<<"############# Check #######################"<<std::endl;
					cs.check();
				}
			} break;
		case 3:/*call CreateSystem::init create, MonteCarlo::run*/
			{
				unsigned int tmax(P.find("tmax",i,false)?P.get<unsigned int>(i):10);
				cs.set_param(&P,NULL);
				std::cout<<"############# Init GenericSystem ##########"<<std::endl;
				cs.construct_GenericSystem(NULL,NULL);
				if(cs.get_status()==2){
					std::cout<<"############# Create GenericSystem ########"<<std::endl;
					cs.create();
					std::cout<<"############# Create MCSystem #############"<<std::endl;
					MCSystem* S(NULL);
					if(cs.use_complex()){
						S = new SystemFermionic<std::complex<double> >(*dynamic_cast<const Fermionic<std::complex<double> >*>(cs.get_system())); 
					} else {
						S = new SystemFermionic<double>(*dynamic_cast<const Fermionic<double>*>(cs.get_system())); 
					}
					std::cout<<"############# Init Monte Carlo ############"<<std::endl;
					MonteCarlo sim(S,tmax);
					sim.thermalize(1e6);
					std::cout<<"############# Run Monte Carlo #############"<<std::endl;
					sim.run();
					S->complete_analysis(1e-5);
					std::cout<<S->get_energy()<<std::endl;
					delete S;
				}
			} break;
		case 4:/*call CreateSystem::(init,create), MonteCarlo::run and save*/
			{
				unsigned int tmax(P.find("tmax",i,false)?P.get<unsigned int>(i):10);
				cs.construct_GenericSystem(NULL,NULL);
				if(cs.get_status()==2){
					std::cout<<"############# Create GenericSystem ########"<<std::endl;
					cs.create();
					std::cout<<"############# Create MCSystem #############"<<std::endl;
					MCSystem* S(NULL);
					if(cs.use_complex()){
						S = new SystemFermionic<std::complex<double> >(*dynamic_cast<const Fermionic<std::complex<double> >*>(cs.get_system())); 
					} else {
						S = new SystemFermionic<double>(*dynamic_cast<const Fermionic<double>*>(cs.get_system())); 
					}
					std::cout<<"############# Init Monte Carlo ############"<<std::endl;
					MonteCarlo sim(S,tmax);
					sim.thermalize(1e6);
					std::cout<<"############# Run Monte Carlo #############"<<std::endl;
					sim.run();
					S->complete_analysis(1e-5);
					std::cout<<S->get_energy()<<std::endl;
					std::cout<<"############# Save MCSystem ###############"<<std::endl;
					IOFiles out("MCSystem.jdbin",true);
					S->write(out);
					delete S;
				}
			} break;
		case 5:/*call CreateSystem::(init,create), MonteCarlo::run and save*/
			{
				unsigned int tmax(P.find("tmax",i,false)?P.get<unsigned int>(i):10);
				MCSystem* S(NULL);
				std::cout<<"############# Load MCSystem ###############"<<std::endl;
				IOFiles in("MCSystem.jdbin",false);
				if(cs.use_complex()){
					S = new SystemFermionic<std::complex<double> >(in); 
				} else {
					S = new SystemFermionic<double>(in); 
				}
				std::cout<<"############# Init Monte Carlo ############"<<std::endl;
				MonteCarlo sim(S,tmax);
				sim.thermalize(10);
				std::cout<<"############# Run Monte Carlo #############"<<std::endl;
				sim.run();
				S->complete_analysis(1e-5);
				std::cout<<S->get_energy()<<std::endl;
				delete S;
			} break;
		case 6:/*call CreateSystem::(init,create), MonteCarlo::run and save*/
			{
				unsigned int tmax(P.find("tmax",i,false)?P.get<unsigned int>(i):10);
				MCSystem* S(NULL);
				std::cout<<"############# Load MCSystem ###############"<<std::endl;
				IOFiles in("MCSystem.jdbin",false);
				if(cs.use_complex()){
					S = new SystemFermionic<std::complex<double> >(in); 
				} else {
					S = new SystemFermionic<double>(in); 
				}
				std::cout<<"############# Init Monte Carlo ############"<<std::endl;
				MonteCarlo sim(S,tmax);
				sim.thermalize(10);
				std::cout<<"############# Run Monte Carlo #############"<<std::endl;
				sim.run();
				S->complete_analysis(1e-5);
				std::cout<<S->get_energy()<<std::endl;
				std::cout<<"############# Save MCSystem ###############"<<std::endl;
				IOFiles out("MCSystem.jdbin",true);
				S->write(out);
				delete S;
			} break;
		default:
			{
				std::cerr<<"check : unknown option 'what', options are :"<<std::endl;
				std::cerr<<"      - init                        : 0"<<std::endl;
				std::cerr<<"      - init + check                : 1"<<std::endl;
				std::cerr<<"      - init + create + check       : 2"<<std::endl;
				std::cerr<<"      - init + create + run         : 3"<<std::endl;
				std::cerr<<"      - init + create + run + write : 4"<<std::endl;
				std::cerr<<"      - load + run                  : 5"<<std::endl;
				std::cerr<<"      - load + run + rewrite        : 6"<<std::endl;
			}
	}
}
