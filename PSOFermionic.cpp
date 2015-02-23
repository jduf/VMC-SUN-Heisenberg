#include "PSOFermionic.hpp"

PSOFermionic::PSOFermionic(Parseur* P):
	PSO(P->get<unsigned int>("Nparticles"),P->get<unsigned int>("Nfreedom"),P->get<double>("cg"),P->get<double>("cp"),P->get<double>("maxiter")),
	tmax_(P->get<unsigned int>("tmax"))
{
	system_.set("N",P->get<unsigned int>("N"));
	system_.set("m",P->get<unsigned int>("m"));
	system_.set("n",P->get<unsigned int>("n"));
	system_.set("bc",P->get<int>("bc"));
	system_.set("type",P->get<unsigned int>("type"));
	system_.set("wf",P->get<std::string>("wf"));
}

double PSOFermionic::f(Vector<double> const& x){
	Container system_param(system_);
	system_param.set("delta",x(0));
	CreateSystem cs(&system_param);
	cs.init();
	if(cs.get_status()==2){
		cs.create();
		if(cs.get_status()==1){
			if( cs.use_complex()){ return monte_carlo<std::complex<double> >(cs,x); } 
			else { return monte_carlo<double>(cs,x); }
		}
	}
	return 0;
}

unsigned int MCSim::cmp_for_fuse(MCSim* list, MCSim* new_elem) { 
	for(unsigned int i(0);i<list->param_.size();i++){
		if(list->param_(i) > new_elem->param_(i)){ return 0; }
		if(list->param_(i) < new_elem->param_(i)){ return 1; }
	}
	return 2;
}

void MCSim::fuse(MCSim* list, MCSim* new_elem) { 
	list->E_.merge(new_elem->E_);
	list->N_++;
}

std::ostream& operator<<(std::ostream& flux, MCSim const& mcsim){
	mcsim.print(flux);
	return flux;
}
