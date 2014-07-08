/*!  @file mc.cpp */

#include "MonteCarlo.hpp"

template<typename Type>
void run(CreateSystem const& cs, unsigned int const& nruns, unsigned int const& tmax);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int nruns(P.get<unsigned int>("nruns"));
	unsigned int tmax(P.get<unsigned int>("tmax"));
	CreateSystem cs(P);
	if(!P.status()){
		do{
			cs.init();
			if(cs.get_status()==1){
				cs.create();
				if(!cs.is_degenerate()){
					if(cs.use_complex()){ run<std::complex<double> >(cs,nruns,tmax); } 
					else { run<double>(cs,nruns,tmax); }
				}
			}
		} while(!cs.is_over());
	}
	return cs.get_status();
}

template<typename Type>
void run(CreateSystem const& cs, unsigned int const& nruns, unsigned int const& tmax){
	IOFiles* file_results(cs.open_and_get_jd_write());
	RST rst;
	rst.title("Input","-");
	file_results->add_to_header(rst.get());
	cs.save(*file_results);

	rst.set();
	rst.title("Simulation's parameters","-");
	file_results->add_to_header(rst.get());
	(*file_results)("number of simulations runned",nruns);
	(*file_results)("tmax",tmax);

	rst.set();
	rst.title("Results","-");
	file_results->add_to_header(rst.get());

	Data<double> E;
	DataSet<double> corr;
	DataSet<double> long_range_corr;
	E.set_conv(true);
	corr.set(cs.get_system()->get_corr().size());
	long_range_corr.set(cs.get_system()->get_long_range_corr().size());

#pragma omp parallel for 
	for(unsigned int i=0;i<nruns;i++){
		Rand rnd(4,omp_get_thread_num());

		MCSystem<Type>* S(NULL);
		if(cs.is_bosonic()){ S = new SystemBosonic<Type>(*dynamic_cast<const Bosonic<Type>*>(cs.get_system()),rnd); } 
		else { S = new SystemFermionic<Type>(*dynamic_cast<const Fermionic<Type>*>(cs.get_system()),rnd); }

		MonteCarlo<Type> sim(S,tmax,rnd);
		sim.run();

#pragma omp critical
		{
			E.add_sample(S->get_energy());
			corr.add_sample(S->get_corr());
			long_range_corr.add_sample(S->get_long_range_corr());
			S->save(*file_results);
		}
		delete S;
	}

	E.complete_analysis();
	corr.complete_analysis();
	long_range_corr.complete_analysis();

	rst.set();
	rst.title("Mean results","-");
	file_results->add_to_header(rst.get());
	(*file_results)("energy per site",E);
	(*file_results)("correlation on links",corr);
	(*file_results)("long range correlation",long_range_corr);
	cs.close_jd_write();
}
