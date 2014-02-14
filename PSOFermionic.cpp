#include "PSOFermionic.hpp"

PSOFermionic::PSOFermionic(Parseur& P,unsigned int Nfreedom):
	PSO(P.get<unsigned int>("Nbees"),Nfreedom,P.get<double>("cg"),P.get<double>("cp"),P.get<double>("maxiter")),
	results_("bla"),
	CS_(P)
{
	std::string wf(P.get<std::string>("wf"));
}

PSOFermionic::~PSOFermionic(){ }

double PSOFermionic::run(Vector<double> const& x){
	Container input;
#pragma omp critical
	{
		CS_.create(x);
		CS_.properties(input);
		CS_.save();
	}

	MonteCarlo<double> sim(input);
	sim.run(1,1e6);
	std::cout<<x<<" "<<sim.get_energy()<<" "<<sim.get_status()<<std::endl;
	return sim.get_energy();
}
