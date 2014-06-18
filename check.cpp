/*!  @file check.cpp */

#include "CreateSystem.hpp"
#include "MonteCarlo.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	CreateSystem cs(P);
	unsigned int i(0);
	unsigned int type(P.get<unsigned int>("type"));
	switch(type){
		case 1:
		case 2:/*run a normal MonteCarlo*/
			{
				unsigned int type(P.get<unsigned int>(i));
				double param(P.get<double>("param"));
				unsigned int tmax(3);
				cs.create(param,type);
				IOFiles w("check.jdbin",true);
				cs.save(w);
				if(cs.use_complex()){
					MCSystem<std::complex<double> >* S(NULL);
					S = new SystemFermionic<std::complex<double> >(cs,type);
					MonteCarlo<std::complex<double> > sim(S,tmax);
					sim.run();
					(sim.get_system())->save(w);
					delete S;
				} else {
					MCSystem<double>* S(NULL);
					S = new SystemFermionic<double>(cs,type);
					MonteCarlo<double> sim(S,tmax);
					sim.run();
					(sim.get_system())->save(w);
					delete S;
				}
			} break;
		case 3:/*call CreateSystem::check*/
			{ cs.check(); } break;
		default:{std::cerr<<"check : unknown type"<<std::endl;}
	}
}
