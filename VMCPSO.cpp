#include "VMCPSO.hpp"

VMCPSO::VMCPSO(Parseur& P):
	Swarm<MCParticle>(P.get<unsigned int>("Nparticles"),P.get<unsigned int>("maxiter"),P.get<unsigned int>("Nfreedom"),P.get<double>("cg"),P.get<double>("cp")),
	VMCMinimization(P,"PSO")
{}

//void VMCPSO::set_x(unsigned int const& i, Vector<double> const& x){
	//VMCMinimization::set_x(i,x);
	//Optimization::set_limit(i,0,x.size());
//}

void VMCPSO::init(bool const& clear_particle_history, bool const& create_particle_history){
	set_time();
	pso_info_.title("New PSO run",'-');

	if(clear_particle_history){ 
		pso_info_.text("clear history"+RST::nl_);
		for(unsigned int p(0);p<Nparticles_;p++){
			std::dynamic_pointer_cast<MCParticle>(particle_[p])->clear_history();
		}
	}

	init_PSO(100); 
	if(clear_particle_history && create_particle_history && all_results_.size()){
		pso_info_.text("set history"+RST::nl_);
		unsigned int size(0);
		while(all_results_.target_next()){
			if(Optimization::within_limit(all_results_.get().get_param())){ size++; }
		}
		unsigned int Npp(size/Nparticles_);
		unsigned int s;
		std::shared_ptr<MCParticle> P;
		for(unsigned int p(0);p<Nparticles_;p++){
			if(p == Nparticles_ - (size%Nparticles_) ){ Npp++; }
			P = std::dynamic_pointer_cast<MCParticle>(particle_[p]);
			s = 0;
			while( s!=Npp && all_results_.target_next()){
				if(Optimization::within_limit(all_results_.get().get_param())){
					s++;
					P->add_to_history(all_results_.get_ptr());
				}
			}
			P->select_new_best();
		}
	}
}

bool VMCPSO::evaluate(unsigned int const& p){
	std::shared_ptr<MCParticle> P(std::dynamic_pointer_cast<MCParticle>(particle_[p]));
	std::shared_ptr<MCSim> sim(compute_vmc(P->get_x()));
	return (sim.get()?P->update(sim):false);
}

void VMCPSO::plot() const {
	std::string filename(get_filename());
	IOFiles data(filename+".dat",true);
	all_results_.set_target();
	while(all_results_.target_next()){
		data<<all_results_.get().get_param()<<all_results_.get().get_S()->get_energy()<<IOFiles::endl;
	}
	all_results_.target_next();
	Gnuplot gp("./",filename);
	unsigned int N(all_results_.get().get_param().size());
	for(unsigned int i(0);i<N;i++){
		gp+=std::string(!i?"plot":"    ")+" '"+filename+".dat' u "+my::tostring(N+1)+":"+my::tostring(i+1)+":"+my::tostring(N+2)+" w xe notitle"+(i==N-1?"":",\\");
	}
	gp.save_file();
	gp.create_image(true);
}

void VMCPSO::print() const {
	Swarm::print();
	std::cout<<"Print whole history ("<< all_results_.size()<<")"<<std::endl;
	while( all_results_.target_next() ){
		std::cout<<all_results_.get_ptr()<<" ";
		all_results_.get().print();
		std::cout<<std::endl;
	}
}

void VMCPSO::save_best(unsigned int const& nsave){
	if(all_results_.size()){
		List<MCSim> best;
		while(all_results_.target_next()){ best.add_sort(all_results_.get_ptr(),VMCPSO::sort_per_energy); }
		best.set_target();
		unsigned int i(0);
		while(best.target_next() && i++<nsave){ best.get().save(&system_param_); }
	} else {
		std::cerr<<"void VMCPSO::save(unsigned int const& nsave) : there is no data"<<std::endl;
	}
}
