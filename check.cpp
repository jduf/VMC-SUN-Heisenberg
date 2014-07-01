/*!  @file check.cpp */

#include "CreateSystem.hpp"
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
				cs.create();
				IOFiles w("check.jdbin",true);
				cs.save(w);
				if(cs.use_complex()){
					MCSystem<std::complex<double> >* S(NULL);
					S = new SystemFermionic<std::complex<double> >(*dynamic_cast<const Fermionic<std::complex<double> >*>(cs.get_system())); 
					MonteCarlo<std::complex<double> > sim(S,tmax);
					sim.run();
					sim.get_system()->save(w);
					delete S;
				} else {
					MCSystem<double>* S(NULL); 
					S = new SystemFermionic<double>(*dynamic_cast<const Fermionic<double>*>(cs.get_system())); 
					MonteCarlo<double> sim(S,tmax);
					sim.run();
					sim.get_system()->save(w);
					delete S;
				}
			} break;
		case 2:/*call CreateSystem::check*/
			{ 
				cs.init();
				cs.check(); 
			} break;
		default:{std::cerr<<"check : unknown what"<<std::endl;}
	}
}
