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
		std::cout<<"ok cs create"<<std::endl;

		IOFiles w("test.jdbin",true);
		cs.save(w);
		if(cs.use_complex()){
			MCSystem<std::complex<double> >* S(NULL);
			S = new SystemFermionic<std::complex<double> >; 
			S->set(cs,type);
			std::cout<<"ok"<<std::endl;
			MonteCarlo<std::complex<double> > sim(S,tmax);
			sim.run();
			(sim.get_system())->save(w);
			delete S;
		} else {
			MCSystem<double>* S(NULL);
			S = new SystemFermionic<double>; 
			S->set(cs,type);
			MonteCarlo<double> sim(S,tmax);
			std::cout<<"ok"<<std::endl;
			sim.run();
			(sim.get_system())->save(w);
			delete S;
		}
	} else {
		//cs.check();
	}
}
