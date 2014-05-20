/*!  @file check.cpp */

#include "CreateSystem.hpp"
#include "MonteCarlo.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	CreateSystem cs(P);
	unsigned int i(0);
	if(P.search("type",i)){
		unsigned int type(P.get<unsigned int>(i));
		double param(P.get<double>("param"));
		unsigned int tmax(3);
		cs.create(param);

		IOFiles w("test.jdbin",true);
		cs.save(w);
		if(cs.use_complex()){
			MCSystem<std::complex<double> >* S(NULL);
			std::cout<<"ok"<<std::endl;
			S = new SystemFermionic<std::complex<double> >(cs.get_system(),cs.get_fermionic<std::complex<double> >()); 
			S->set_type(type);
			MonteCarlo<std::complex<double> > sim(S,tmax);
			sim.run();
			(sim.get_system())->save(w);
			delete S;
		} else {
			MCSystem<double>* S(NULL);
			std::cout<<"ok"<<std::endl;
			S = new SystemFermionic<double>(cs.get_system(),cs.get_fermionic<double>()); 
			S->set_type(type);
			MonteCarlo<double> sim(S,tmax);
			sim.run();
			(sim.get_system())->save(w);
			delete S;
		}
	} else {
		//cs.check();
	}
}
