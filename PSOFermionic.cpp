#include "PSOFermionic.hpp"

PSOFermionic::PSOFermionic(Parseur& P,unsigned int Nfreedom, unsigned int N_MC):
	PSO(P.get<unsigned int>("Nbees"),Nfreedom,P.get<double>("cg"),P.get<double>("cp"),P.get<double>("maxiter")),
	CS_(P),
	N_MC_(N_MC)
{ 
}

PSOFermionic::~PSOFermionic(){}

double PSOFermionic::run(Vector<double> const& x){
	Vector<double> t(x.size()+1);
	t(0) = 1.0;
	for(unsigned int i(1); i<t.size();i++){
		t(i) = x(i-1);
	}
	std::cerr<<"PSOFermionic is wrong"<<std::endl;
	CreateSystem CS(CS_,1);

	MonteCarlo<double> sim(CS,N_MC_);
	sim.run();
	return sim.get_energy();
}
