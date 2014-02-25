#include "PSOFermionic.hpp"

PSOFermionic::PSOFermionic(Parseur& P,unsigned int Nfreedom, unsigned int N_MC):
	PSO(P.get<unsigned int>("Nbees"),Nfreedom,P.get<double>("cg"),P.get<double>("cp"),P.get<double>("maxiter")),
	N_MC_(N_MC)
{ 
	CreateSystem(P,param_);
}

PSOFermionic::~PSOFermionic(){}

double PSOFermionic::run(Vector<double> const& x){
	Container input;
	CreateSystem CS(param_);
	Vector<double> t(x.size()+1);
	t(0) = 1.0;
	for(unsigned int i(1); i<t.size();i++){
		t(i) = x(i-1);
	}
	CS.create(t);
	CS.properties(input);

	MonteCarlo<double> sim(input);
	sim.run(2,N_MC_);
	std::cout<<x<<" "<<sim.get_energy()<<" "<<sim.get_status()<<std::endl;
	return sim.get_energy();
}
