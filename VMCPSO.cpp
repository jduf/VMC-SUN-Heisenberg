#include "VMCPSO.hpp"

VMCPSO::VMCPSO(Parseur* P):
	Swarm<MCParticle>(P->get<unsigned int>("Nparticles"),P->get<unsigned int>("maxiter"),P->get<unsigned int>("Nfreedom"),P->get<double>("cg"),P->get<double>("cp")),
	tmax_(P->get<unsigned int>("tmax"))
{
	system_param_.set("N",P->get<unsigned int>("N"));
	system_param_.set("m",P->get<unsigned int>("m"));
	system_param_.set("n",P->get<unsigned int>("n"));
	system_param_.set("bc",P->get<int>("bc"));
	system_param_.set("wf",P->get<std::string>("wf"));
}

bool VMCPSO::evaluate(unsigned int const& p){
	std::shared_ptr<MCParticle> P(std::dynamic_pointer_cast<MCParticle>(particle_[p]));
	std::shared_ptr<MCSim> sim(std::make_shared<MCSim>(P->get_x()));
	bool tmp_test;
#pragma omp critical(all_results_)
	{
		if(all_results_.find_sorted(sim,MCSim::cmp_for_fuse)){ 
			tmp_test = true;
			sim->copy_S(all_results_.get().get_S()); 
			sim->get_S()->set(); //to clear the binning
		} else {
			tmp_test = false;
			sim->create_S(&system_param_,&P->get_x());
		}
	}
	if(sim->is_created()){
		MonteCarlo mc(sim->get_S().get(),tmax_);
		mc.thermalize(tmp_test?10:1e6);
		mc.run();

#pragma omp critical(all_results_)
		{
			if(all_results_.find_sorted(sim,MCSim::cmp_for_fuse)){ 
				all_results_.fuse_with_target(sim,MCSim::fuse); 
				tmp_test = P->update(all_results_.get_ptr()); 
			} else {
				all_results_.add_after_target(sim); 
				tmp_test = P->update(sim); 
			}
		}
		return tmp_test;
	} else {
		std::cerr<<"bool VMCPSO::evaluate(unsigned int const& p) : not valid parameter : "<<particle_[p]->get_x()<<std::endl;
		return false;
	}
}

void VMCPSO::refine(unsigned int const& Nrefine, double const& tol, unsigned int const& tmax){
	auto sort = [](MCSim const& a, MCSim const& b){ 
		return a.get_S()->get_energy().get_x()<b.get_S()->get_energy().get_x();
	};

	List<MCSim> best;
	all_results_.set_target();
	while(all_results_.target_next()){
		best.add_sort(all_results_.get_ptr(),sort);
	}
	unsigned int N(all_results_.size());
	N  = (N<Nrefine?N:Nrefine);
	best.set_target();
#pragma omp parallel for schedule(dynamic,1)
	for(unsigned int i=0;i<N;i++){
		std::shared_ptr<MCSim> sim;
#pragma omp critical
		{
			best.target_next();
			sim = best.get_ptr();
		}
		MonteCarlo mc(sim->get_S().get(),tmax);
		while(sim->get_S()->get_energy().get_dx()>tol){
			mc.run();
			sim->get_S()->complete_analysis(1e-5);
		}
#pragma omp critical
		{
			std::cout<<i<<" sim refined "<<sim->get_S()->get_energy()<<std::endl;
		}
	}
	//best.set_free();
	//while(best.target_next()){
	//std::cout<<best.get().get_S()->get_energy()<<std::endl;
	//}
}

void VMCPSO::plot() const {
	IOFiles data("data.dat",true);
	all_results_.set_target();
	while(all_results_.target_next()){
		data<<all_results_.get().get_param()<<all_results_.get().get_S()->get_energy()<<IOFiles::endl;
	}
	all_results_.target_next();
	Gnuplot gp("./","test");
	unsigned int N(all_results_.get().get_param().size());
	gp+="plot 'data.dat' u "+my::tostring(N+1)+":1:"+my::tostring(N+2)+" w xe"+(N==1?"":",\\"); 
	for(unsigned int i(1);i<N;i++){
		gp+="     'data.dat' u "+my::tostring(N+1)+":"+my::tostring(i+1)+":"+my::tostring(N+2)+" w xe"+(i==N-1?"":",\\");
	}
	gp.save_file();
	gp.create_image(true);
}

void VMCPSO::print() const {
	Swarm::print();
	std::cout<<"Print whole history"<<std::endl;
	all_results_.set_target();
	while( all_results_.target_next() ){
		std::cout<<all_results_.get_ptr()<<" ";
		all_results_.get().print();
		std::cout<<std::endl;
	}
}

void VMCPSO::write(IOFiles& w) const {
	w<<all_results_.size();
	while(all_results_.target_next()){ all_results_.get().write(w); }
}

void VMCPSO::read(IOFiles& r, bool create_particle_history){
	unsigned int size;
	r>>size;
	while(size--){ all_results_.add_end(std::make_shared<MCSim>(r)); }
	//if(create_particle_history){
		//size = 0;
		//while(all_results_.target_next()){
			//if(Optimization::within_limit(all_results_.get().get_param())){ size++; }
		//}
		//unsigned int Npp(size/Nparticles_);
		//unsigned int s;
		//std::shared_ptr<MCParticle> P;
		//for(unsigned int p(0);p<Nparticles_;p++){
			//if(p==Nparticles_-size%Nparticles_){ Npp++; }
			//P = std::dynamic_pointer_cast<MCParticle>(particle_[p]);
			//s = 0;
			//while( s!=Npp && all_results_.target_next()){
				//if(Optimization::within_limit(all_results_.get().get_param())){
					//s++;
					//P->add_to_history(all_results_.get_ptr());
				//}
			//}
			//P->select_new_best();
		//}
	//}
}

void VMCPSO::complete_analysis(double const& tol){
	all_results_.set_target();
	while ( all_results_.target_next() ){
		all_results_.get().get_S()->complete_analysis(tol);
	}
}
