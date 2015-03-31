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
	list_elem.E_.merge(new_elem.E_);
	list_elem.N_++;
	//new_elem.E_ = list_elem.E_;
	new_elem = list_elem;
}

std::ostream& operator<<(std::ostream& flux, MCSim const& mcsim){
	mcsim.print(flux);
	return flux;
}
/*}*/

/*{MCParticle*/
unsigned int MCParticle::pos_iter = 0;
Vector<double>MCParticle::pos = Vector<double>(8);

void MCParticle::move(Vector<double> const& bx_all){
	//Particle::move(bx_all);
	//unsigned int n;
	//double dx(0.01);
	//for(unsigned int j(0);j<Nfreedom_;j++){
		//n=0;
		//if(std::abs(x_(j))<dx/2){ n=1; }
		//if(std::abs(x_(j)-min_(j))<dx/2){ n=2; }
		//if(std::abs(x_(j)-max_(j))<dx/2){ n=3; }
		//switch(n){
			//case 0:{ x_(j) = std::round(x_(j)/dx)*dx; }break;
			//case 1:{ x_(j) = 0; }break;
			//case 2:{ x_(j) = min_(j); }break;
			//case 3:{ x_(j) = max_(j); }break;
		//}
	//}
	
	(void)(bx_all);
#pragma omp critical(update_pos_iter)
	{
		x_(0) = MCParticle::pos(MCParticle::pos_iter);
		std::cout<<"pos_iter="<<MCParticle::pos_iter<<" x_="<<x_(0)<<std::endl;
		MCParticle::pos_iter++; 
	}
}

void MCParticle::update_particle_history(std::shared_ptr<MCSim>& new_elem){
	if(history_.find(new_elem,MCSim::cmp_for_fuse)){ history_.add_after_move(new_elem); }
}

void MCParticle::print(std::ostream& flux) const {
	Particle::print(flux);
	flux<<std::endl<<"particle history"<<std::endl;
	do{
		flux<<history_.get_ptr()<<" "<<history_.get().get_param()<<" "<<history_.get().get_energy()<<" "<<history_.get().get_N()<<IOFiles::endl;
	} while ( history_.move_forward() );
}
/*}*/

/*{PSOFermionic*/
PSOFermionic::PSOFermionic(Parseur* P):
	Swarm<MCParticle>(P->get<unsigned int>("Nparticles"),P->get<unsigned int>("maxiter"),P->get<unsigned int>("Nfreedom"),P->get<double>("cg"),P->get<double>("cp")),
	tmax_(P->get<unsigned int>("tmax"))
{
	MCParticle::pos(0) = 0.1;
	MCParticle::pos(1) = 0.2;
	MCParticle::pos(2) = 0.3;
	MCParticle::pos(3) = 0.1;
	MCParticle::pos(4) = 0.2;
	MCParticle::pos(5) = 0.2;
	MCParticle::pos(6) = 0.1;
	MCParticle::pos(7) = 0.4;
	system_.set("N",P->get<unsigned int>("N"));
	system_.set("m",P->get<unsigned int>("m"));
	system_.set("n",P->get<unsigned int>("n"));
	system_.set("bc",P->get<int>("bc"));
	system_.set("wf",P->get<std::string>("wf"));
}

bool PSOFermionic::is_better_x(unsigned int const& p){
	Container system_param(system_);
	//Vector<double> t(10,1);
	//for(unsigned int i(0);i<t.size()-1;i++){ t(i) = x(i); }
	//Vector<double> mu(5,0);
	//for(unsigned int i(0);i<mu.size();i++){ mu(i) = x(i+t.size()-1); }
	//system_param.set("t",t);
	//system_param.set("mu",mu);
	//system_param.set("td",x(0));
	//
	system_param.set("delta",p_[p]->get_x()(0));

	//Vector<double> mu(1,0);
	//system_param.set("mu",mu);
	//Vector<double> phi(1,M_PI/4.0);
	//system_param.set("phi",phi);
	//system_param.set("t",p_[p]->get_x());

	//system_param.set("nu",x);

	CreateSystem cs(&system_param);
	cs.init();
	if(cs.get_status()==2){
		cs.create();
		if(cs.get_status()==1){
			if( cs.use_complex()){ return monte_carlo<std::complex<double> >(cs,p); } 
			else { return monte_carlo<double>(cs,p); }
		}
	}
	return false;
}

void PSOFermionic::plot(){
	IOFiles data("data.dat",true);
	do{
		data<<all_results_.get().get_param()<<" "<<all_results_.get().get_energy()<<" "<<all_results_.get().get_N()<<IOFiles::endl;
	} while ( all_results_.move_forward() );
	Gnuplot gp("./","test");
	gp+="plot 'data.dat' u 1:2:3 w e";
	gp.save_file();
	gp.create_image(true);
}

void PSOFermionic::print(std::ostream& flux) const {
	Swarm::print(flux);
	std::cout<<"Print whole history"<<std::endl;
	do{
		flux<<all_results_.get_ptr()<<" "<<all_results_.get().get_param()<<" "<<all_results_.get().get_energy()<<" "<<all_results_.get().get_N()<<IOFiles::endl;
	} while ( all_results_.move_forward() );
}

std::ostream& operator<<(std::ostream& flux, PSOFermionic const& pso){
	pso.print(flux);
	return flux;
}
/*}*/

