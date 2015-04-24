#include "PSOMonteCarlo.hpp"

PSOMonteCarlo::PSOMonteCarlo(Parseur* P):
	Swarm<MCParticle>(P->get<unsigned int>("Nparticles"),P->get<unsigned int>("maxiter"),P->get<unsigned int>("Nfreedom"),P->get<double>("cg"),P->get<double>("cp")),
	tmax_(P->get<unsigned int>("tmax"))
{
	system_.set("N",P->get<unsigned int>("N"));
	system_.set("m",P->get<unsigned int>("m"));
	system_.set("n",P->get<unsigned int>("n"));
	system_.set("bc",P->get<int>("bc"));
	system_.set("wf",P->get<std::string>("wf"));
}

bool PSOMonteCarlo::is_better_x(unsigned int const& p){
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
			sim->create_S(&system_,&p_[p]->get_x());
		}
	}
	if(sim->is_created()){
		MonteCarlo mc(sim->get_S().get(),tmax_);
		mc.thermalize(tmp_test?10:1e6);
		mc.run();

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
	} else {
		std::cerr<<"bool PSOMonteCarlo::is_better_x(unsigned int const& p) : not valid parameter : "<<p_[p]->get_x()<<std::endl;
		return false;
	}
}

void PSOMonteCarlo::plot() const {
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

void PSOMonteCarlo::print() const {
	Swarm::print();
	std::cout<<"Print whole history"<<std::endl;
	all_results_.set_free();
	while( all_results_.go_to_next() ){
		std::cout<<all_results_.get_ptr()<<" ";
		all_results_.get().print();
	}
}

void PSOMonteCarlo::write(IOFiles& w) const {
	w<<all_results_.size();
	while(all_results_.go_to_next()){ all_results_.get().write(w); }
}

void PSOMonteCarlo::read(IOFiles& r){
	unsigned int size;
	r>>size;
	while(size--){ all_results_.add_end(std::make_shared<MCSim>(r)); }
}

void PSOMonteCarlo::complete_analysis(double tol){
	all_results_.set_free();
	while ( all_results_.go_to_next() ){
		all_results_.get().get_S()->complete_analysis(tol);
	}

}
