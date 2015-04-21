#include "PSOFermionic.hpp"

/*{MCSim*/
unsigned int MCSim::cmp_for_fuse(MCSim const& list, MCSim const& new_elem) {
	for(unsigned int i(0);i<list.param_.size();i++){
		if(list.param_(i) > new_elem.param_(i)){ return 0; }
		if(list.param_(i) < new_elem.param_(i)){ return 1; }
	}
	return 2;
}

void MCSim::fuse(MCSim& list_elem, MCSim& new_elem) { 
	list_elem.get_S().get()->get_energy().merge(new_elem.get_S().get()->get_energy());
}

void MCSim::create_S(Container* C){
	CreateSystem cs(C);
	cs.init();
	if(cs.get_status()==2){
		cs.create();
		if(cs.get_status()==1){
			if( cs.use_complex()){
				if(cs.is_bosonic()){
					S_.reset(new SystemBosonic<std::complex<double> >(*dynamic_cast<const Bosonic<std::complex<double> >*>(cs.get_system())));
				} else {
					S_.reset(new SystemFermionic<std::complex<double> >(*dynamic_cast<const Fermionic<std::complex<double> >*>(cs.get_system())));
				}
			} else {
				if(cs.is_bosonic()){
					S_.reset(new SystemBosonic<double>(*dynamic_cast<const Bosonic<double>*>(cs.get_system())));
				} else {
					S_.reset(new SystemFermionic<double>(*dynamic_cast<const Fermionic<double>*>(cs.get_system())));
				}
			}
		}
	}
}

void MCSim::copy_S(std::unique_ptr<MCSystem> const& S){
	S_ = S->clone();
}

std::ostream& operator<<(std::ostream& flux, MCSim const& mcsim){
	mcsim.print(flux);
	return flux;
}
/*}*/

/*{MCParticle*/
//unsigned int MCParticle::pos_iter = 0;
//Vector<double>MCParticle::pos = Vector<double>(16);

void MCParticle::move(Vector<double> const& bx_all){
	Particle::move(bx_all);
	unsigned int n;
	double dx(0.01);
	for(unsigned int j(0);j<Nfreedom_;j++){
		n=0;
		if(std::abs(x_(j))<dx/2){ n=1; }
		if(std::abs(x_(j)-min_(j))<dx/2){ n=2; }
		if(std::abs(x_(j)-max_(j))<dx/2){ n=3; }
		switch(n){
			case 0:{ x_(j) = std::round(x_(j)/dx)*dx; }break;
			case 1:{ x_(j) = 0; }break;
			case 2:{ x_(j) = min_(j); }break;
			case 3:{ x_(j) = max_(j); }break;
		}
	}

	//(void)(bx_all);
	//#pragma omp critical(update_pos_iter)
	//{
	//x_(0) = MCParticle::pos(MCParticle::pos_iter);
	//MCParticle::pos_iter++; 
	//}
}

bool MCParticle::update(std::shared_ptr<MCSim> const& new_elem){
	if(history_.find_sorted(new_elem,MCSim::cmp_for_fuse)){ history_.set_free(); }
	else{ history_.add_after_free(new_elem); }

	double tmp;
	bool is_better(false);
	while(history_.go_to_next()){
		tmp = history_.get().get_S()->get_energy().get_x();
		if(tmp<fbx_){
			bx_ = history_.get().get_param();
			fbx_ = tmp;
			is_better = true;
		}
	}
	return is_better;
}

void MCParticle::print(std::ostream& flux) const {
	Particle::print(flux);
	flux<<std::endl<<"particle history "<<std::endl;
	history_.set_free();
	while( history_.go_to_next() ){
		flux<<history_.get_ptr()<<" "<<history_.get()<<IOFiles::endl;
	}
}
/*}*/

/*{PSOFermionic*/
PSOFermionic::PSOFermionic(Parseur* P):
	Swarm<MCParticle>(P->get<unsigned int>("Nparticles"),P->get<unsigned int>("maxiter"),P->get<unsigned int>("Nfreedom"),P->get<double>("cg"),P->get<double>("cp")),
	tmax_(P->get<unsigned int>("tmax"))
{
	//MCParticle::pos(0) = 0.1;
	//MCParticle::pos(1) = 0.2;
	//MCParticle::pos(2) = 0.3;
	//MCParticle::pos(3) = 0.1;
	//MCParticle::pos(4) = 0.2;
	//MCParticle::pos(5) = 0.2;
	//MCParticle::pos(6) = 0.1;
	//MCParticle::pos(7) = 0.4;
	//MCParticle::pos(8) = 0.6;
	//MCParticle::pos(9) = 0.5;
	//MCParticle::pos(10) = 0.5;
	//MCParticle::pos(11) = 0.5;
	//MCParticle::pos(12) = 0.5;
	//MCParticle::pos(13) = 0.2;
	//MCParticle::pos(14) = 0.1;
	//MCParticle::pos(15) = 0.3;
	system_.set("N",P->get<unsigned int>("N"));
	system_.set("m",P->get<unsigned int>("m"));
	system_.set("n",P->get<unsigned int>("n"));
	system_.set("bc",P->get<int>("bc"));
	system_.set("wf",P->get<std::string>("wf"));
}

bool PSOFermionic::is_better_x(unsigned int const& p){
	std::shared_ptr<MCSim> sim(std::make_shared<MCSim>(p_[p]->get_x()));
	bool tmp_test;
#pragma omp critical(all_results_)
	{
		if(all_results_.find_sorted(sim,MCSim::cmp_for_fuse)){ 
			tmp_test = true;
			sim->copy_S(all_results_.get().get_S()); 
			sim->get_S()->set(); //to clear the binning
		} else {
			tmp_test = false;
			Container system_param(system_);
			system_param.set("t2",p_[p]->get_x()(0));
			sim->create_S(&system_param);
		}
	}

	MonteCarlo mc(sim->get_S().get(),tmax_);
	mc.thermalize(tmp_test?10:1e6);
	mc.run();
	//mc.complete_analysis(1e-5);

	std::shared_ptr<MCParticle> P(std::dynamic_pointer_cast<MCParticle>(p_[p]));
#pragma omp critical(all_results_)
	{
		if(all_results_.find_sorted(sim,MCSim::cmp_for_fuse)){ 
			all_results_.fuse_with_free(sim,MCSim::fuse); 
			tmp_test = P->update(all_results_.get_ptr()); 
		} else {
			all_results_.add_after_free(sim); 
			tmp_test = P->update(sim); 
		}
	}
	return tmp_test;
}

void PSOFermionic::plot() const {
	IOFiles data("data.dat",true);
	all_results_.set_free();
	while(all_results_.go_to_next()){
		data<<all_results_.get().get_param()<<all_results_.get().get_S()->get_energy()<<IOFiles::endl;
	}
	Gnuplot gp("./","test");
	gp+="plot 'data.dat' u 1:2:3 w e";
	gp.save_file();
	gp.create_image(true);
}

void PSOFermionic::print(std::ostream& flux) const {
	Swarm::print(flux);
	std::cout<<"Print whole history"<<std::endl;
	all_results_.set_free();
	while ( all_results_.go_to_next() ){
		flux<<all_results_.get_ptr()<<" "<<all_results_.get()<<IOFiles::endl;
	}
}

void PSOFermionic::complete_analysis(double tol){
	all_results_.set_free();
	while ( all_results_.go_to_next() ){
		all_results_.get().get_S()->complete_analysis(tol);
	}

}

std::ostream& operator<<(std::ostream& flux, PSOFermionic const& pso){
	pso.print(flux);
	return flux;
}
/*}*/
