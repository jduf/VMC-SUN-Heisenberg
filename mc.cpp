/*!  @file mc.cpp */

#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"

template<typename Type>
void run(CreateSystem const& cs, unsigned int const& nruns, unsigned int const& tmax);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int nruns(P.get<unsigned int>("nruns"));
	unsigned int tmax(P.get<unsigned int>("tmax"));
	System s(P);
	CreateSystem cs(&s);
	cs.init(NULL,&P);
	if(!P.locked()){
		if(cs.get_status()==2){
			cs.create();
			if(cs.get_status()==1){
				if(cs.use_complex()){ run<std::complex<double> >(cs,nruns,tmax); }
				else { run<double>(cs,nruns,tmax); }
			}
		}
	}
	return cs.get_status();
}

template<typename Type>
void run(CreateSystem const& cs, unsigned int const& nruns, unsigned int const& tmax){
	(void)(nruns);
	Linux command;
	command.mkdir(cs.get_path());
	IOFiles out(cs.get_path() + cs.get_filename()+".jdbin",true);

	cs.save_param(out);
	cs.get_system()->save_input(out);


	MCSystem* S(NULL);
	if(cs.is_bosonic())
	{ S = new SystemBosonic<Type>
		(*dynamic_cast<const Bosonic<Type>*>(cs.get_system())); } 
	else 
	{ S = new SystemFermionic<Type>
		(*dynamic_cast<const Fermionic<Type>*>(cs.get_system())); }

	S->set_observable(2);
	MonteCarlo sim(S,tmax);
	sim.thermalize(1e6);
	sim.run();
	S->complete_analysis(1e-5);
	S->delete_binning();
	S->save_output(out);
	std::cout<<S->get_energy()<<std::endl;
}
