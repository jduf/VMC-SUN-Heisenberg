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
	cs.set_param(&P,NULL);
	if(!P.locked()){
		cs.construct_GenericSystem(NULL,NULL);
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
	Linux command;
	command.mkdir(cs.get_path());
	IOFiles file_results(cs.get_path() + cs.get_filename()+".jdbin",true);
	cs.get_system()->save_input(file_results);
	cs.save_param(file_results);

	RST rst;
	rst.title("Simulation's parameters",'-');
	file_results.add_header()->add(rst.get());
	file_results.write("number of simulations runned",nruns);
	file_results.write("tmax",tmax);
	std::cout<<file_results.get_header()<<std::endl;
	rst.set();
	rst.title("Results",'-');
	file_results.add_header()->add(rst.get());

	Data<double> E;
	DataSet<double> corr;
	DataSet<double> lr_corr;
	corr.set(cs.get_system()->get_corr().size());
	lr_corr.set(cs.get_system()->get_lr_corr().size());

#pragma omp parallel for 
	for(unsigned int i=0;i<nruns;i++){
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
		std::cout<<S->get_energy()<<std::endl;

#pragma omp critical
		{
			E.add_sample(S->get_energy());
			corr.add_sample(S->get_corr());
			lr_corr.add_sample(S->get_lr_corr());
			file_results.write("energy per site",S->get_energy());
			file_results.write("correlation on links",S->get_corr());
			file_results.write("long range correlation",S->get_lr_corr());
		}
		delete S;
	}

	E.complete_analysis();
	corr.complete_analysis();
	lr_corr.complete_analysis();

	rst.set();
	rst.title("Mean results",'-');
	file_results.add_header()->add(rst.get());
	file_results.write("energy per site",E);
	file_results.write("correlation on links",corr);
	file_results.write("long range correlation",lr_corr);
	std::cout<<E<<std::endl;
}
