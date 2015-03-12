/*!  @file check.cpp */

#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	CreateSystem cs(&P);
	switch(P.get<unsigned int>("what")){
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
		case 3:/*call CreateSystem::init create and check, then run MonteCarlo*/
			{
				unsigned int tmax(10);
				cs.init();
				if(cs.get_status()==2){
					cs.create();
					IOFiles w("check.jdbin",true);
					cs.init_output_file(w);
					cs.save();
					if(cs.use_complex()){
						MCSystem<std::complex<double> >* S(NULL);
						S = new SystemFermionic<std::complex<double> >(*dynamic_cast<const Fermionic<std::complex<double> >*>(cs.get_system())); 
						MonteCarlo<std::complex<double> > sim(S,tmax);
						sim.run();
						w.write("energy per site",S->get_energy());
						w.write("correlation on links",S->get_corr());
						w.write("long range correlation",S->get_lr_corr());
						std::cout<<S->get_energy()<<std::endl;
						delete S;
					} else {
						MCSystem<double>* S(NULL); 
						S = new SystemFermionic<double>(*dynamic_cast<const Fermionic<double>*>(cs.get_system())); 
						MonteCarlo<double> sim(S,tmax);
						sim.run();
						w.write("energy per site",S->get_energy());
						w.write("correlation on links",S->get_corr());
						w.write("long range correlation",S->get_lr_corr());
						std::cout<<S->get_energy()<<std::endl;
						delete S;
					}
				}
			} break;
		default:{std::cerr<<"check : unknown what"<<std::endl;}
	}
}
