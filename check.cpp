/*!@file check.cpp */

#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	std::cout<<"############# Init System #################"<<std::endl;
	System s(P);
	std::cout<<"############# Init CreateSystem ###########"<<std::endl;
	CreateSystem cs(&s);
	unsigned int i(0);
	switch(P.find("what",i,true)?P.get<unsigned int>(i):666){
		case 0:/*call CreateSystem::init*/
			{ 
				std::cout<<"############# Init GenericSystem ##########"<<std::endl;
				cs.init(NULL,&P);
			} break;
		case 1:/*call CreateSystem::init and check*/
			{ 
				std::cout<<"############# Init GenericSystem ##########"<<std::endl;
				cs.init(NULL,&P);
				std::cout<<"############# Check #######################"<<std::endl;
				if(cs.get_status()==2){ cs.check(); }
			} break;
		case 2:/*call CreateSystem::init create and check*/
			{ 
				std::cout<<"############# Init GenericSystem ##########"<<std::endl;
				cs.init(NULL,&P);
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
				std::cout<<"############# Init GenericSystem ##########"<<std::endl;
				cs.init(NULL,&P);
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
					S->set_observable(1);
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
				std::cout<<"############# Init GenericSystem ##########"<<std::endl;
				cs.init(NULL,&P);
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
		case 7:/*check symmetries*/
			{
				Vector<double> t_ref(P.get<std::vector<double> >("t"));
				cs.init(NULL,&P);
				std::vector<Matrix<int> > all_sym;
				cs.get_wf_symmetries(all_sym);
				for(unsigned int j(0);j<all_sym.size();j++){
					Vector<double> t(t_ref);
					Matrix<int> sym(all_sym[j]);
					for(unsigned int i(0);i<sym.row();i++){
						if(sym(i,1)<0){
							t(sym(i,0)) = sym(i,2)*0.1;
						} else {
							t(sym(i,0)) = sym(i,2)*t(sym(i,1));
						}
					}
					//std::cout<<sym<<std::endl<<std::endl;
					std::cout<<"sim["<<j<<"]"<<" "<<i<<" -> "<<t<<std::endl;
				}
			} break;
		default:
			{
				std::cerr<<__PRETTY_FUNCTION__<<" : unknown option 'what', options are :"<<std::endl;
				std::cerr<<"    - init                        : 0"<<std::endl;
				std::cerr<<"    - init + check                : 1"<<std::endl;
				std::cerr<<"    - init + create + check       : 2"<<std::endl;
				std::cerr<<"    - init + create + run         : 3"<<std::endl;
				std::cerr<<"    - init + create + run + write : 4"<<std::endl;
				std::cerr<<"    - load + run                  : 5"<<std::endl;
				std::cerr<<"    - load + run + rewrite        : 6"<<std::endl;
				std::cerr<<"    - load + check_symmetries     : 7"<<std::endl;
			}
	}
}
