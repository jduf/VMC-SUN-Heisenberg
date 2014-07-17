/*!  @file check.cpp */

#include "MonteCarlo.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	CreateSystem cs(P);
	unsigned int what(P.get<unsigned int>("what"));
	switch(what){
		case 1:/*run a normal MonteCarlo*/
			{
				unsigned int tmax(3);
				cs.init();
				if(cs.get_status()==1){
					cs.create();
					IOFiles w("check.jdbin",true);
					cs.init_output_file(w);
					cs.save();
					Rand rnd(4);
					if(cs.use_complex()){
						MCSystem<std::complex<double> >* S(NULL);
						S = new SystemFermionic<std::complex<double> >(*dynamic_cast<const Fermionic<std::complex<double> >*>(cs.get_system()),rnd); 
						MonteCarlo<std::complex<double> > sim(S,tmax,rnd);
						sim.run();
						w("energy per site",S->get_energy());
						w("correlation on links",S->get_corr());
						w("long range correlation",S->get_long_range_corr());
						std::cout<<S->get_energy()<<std::endl;
						delete S;
					} else {
						MCSystem<double>* S(NULL); 
						S = new SystemFermionic<double>(*dynamic_cast<const Fermionic<double>*>(cs.get_system()),rnd); 
						MonteCarlo<double> sim(S,tmax,rnd);
						sim.run();
						w("energy per site",S->get_energy());
						w("correlation on links",S->get_corr());
						w("long range correlation",S->get_long_range_corr());
						std::cout<<S->get_energy()<<std::endl;
						delete S;
					}
				}
			} break;
		case 2:/*call CreateSystem::check*/
			{ 
				cs.init();
				if(cs.get_status()==1){ cs.check(); }
			} break;
		default:{std::cerr<<"check : unknown what"<<std::endl;}
	}
}
