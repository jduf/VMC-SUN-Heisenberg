#include "VMCPSO.hpp"

VMCPSO::VMCPSO(Parseur& P, Minimization& m):
	Swarm<MCParticle>(P.get<unsigned int>("Nparticles"),P.get<unsigned int>("maxiter"),P.get<unsigned int>("Nfreedom"),P.get<double>("cg"),P.get<double>("cp")),
	VMCMinimization(m,"PSO")
{}

void VMCPSO::init(bool const& clear_particle_history, bool const& create_particle_history){
	set_time();
	m_.pso_info_.title("New PSO run",'-');

	std::cout<<"#######################"<<std::endl;
	std::cout<<"#new VMCPSO"<<std::endl;
	std::cout<<"#"<<get_filename()<<std::endl;
	std::cout<<"#contains "<<m_.all_results_.size()<<" samples"<<std::endl;

	if(clear_particle_history){ 
		m_.pso_info_.text("clear history"+RST::nl_);
		for(unsigned int p(0);p<Nparticles_;p++){
			std::dynamic_pointer_cast<MCParticle>(particle_[p])->clear_history();
		}
	}

	std::shared_ptr<MCParticle> MCP;
	for(unsigned int i(0);i<Nparticles_;i++){
		MCP = std::dynamic_pointer_cast<MCParticle>(particle_[i]);
		MCP->set_ps(m_.ps_);
	}
	init_PSO(100); 
	if(clear_particle_history && create_particle_history && m_.all_results_.size()){
		m_.pso_info_.text("set history"+RST::nl_);
		unsigned int size(0);
		while(m_.all_results_.target_next()){
			if(m_.within_limit(m_.all_results_.get().get_param())){ size++; }
		}
		unsigned int Npp(size/Nparticles_);
		unsigned int s;
		std::shared_ptr<MCParticle> MCP;
		for(unsigned int p(0);p<Nparticles_;p++){
			if(p == Nparticles_ - (size%Nparticles_) ){ Npp++; }
			MCP = std::dynamic_pointer_cast<MCParticle>(particle_[p]);
			s = 0;
			while( s!=Npp && m_.all_results_.target_next()){
				if(m_.within_limit(m_.all_results_.get().get_param())){
					MCP->add_to_history(m_.all_results_.get_ptr());
					s++;
				}
			}
			MCP->select_new_best();
		}
	}
}

void VMCPSO::set_ps(unsigned int const& i, Vector<double> const& ps){
	m_.set_ps(i,ps);
	if(i<m_.Nfreedom_){
		Particle::set_limit(i,0,m_.ps_[i].size());
	} else {
		std::cerr<<"void VMCPSO::set_x(unsigned int const& i, Vector<double> const& x) : i>=Nfreedom"<<std::endl;
	}
}

bool VMCPSO::evaluate(unsigned int const& p){
	std::shared_ptr<MCParticle> MCP(std::dynamic_pointer_cast<MCParticle>(particle_[p]));
	std::shared_ptr<MCSim> sim(compute_vmc(MCP->get_param()));
	return (sim.get()?MCP->update(sim):false);
}

void VMCPSO::plot() const {
	std::string filename(get_filename());
	IOFiles data(filename+".dat",true);
	m_.all_results_.set_target();
	while(m_.all_results_.target_next()){
		data<<m_.all_results_.get().get_param()<<m_.all_results_.get().get_S()->get_energy()<<IOFiles::endl;
	}
	m_.all_results_.target_next();
	Gnuplot gp("./",filename);
	unsigned int N(m_.all_results_.get().get_param().size());
	for(unsigned int i(0);i<N;i++){
		gp+=std::string(!i?"plot":"    ")+" '"+filename+".dat' u "+my::tostring(N+1)+":"+my::tostring(i+1)+":"+my::tostring(N+2)+" w xe notitle"+(i==N-1?"":",\\");
	}
	gp.save_file();
	gp.create_image(true);
}

void VMCPSO::print() const {
	Swarm::print();
	std::cout<<"Print whole history ("<< m_.all_results_.size()<<")"<<std::endl;
	while( m_.all_results_.target_next() ){
		std::cout<<m_.all_results_.get_ptr()<<" ";
		m_.all_results_.get().print();
		std::cout<<std::endl;
	}
}

void VMCPSO::save_best(unsigned int const& nsave){
	if(m_.all_results_.size()){
		List<MCSim> best;
		while(m_.all_results_.target_next()){ best.add_sort(m_.all_results_.get_ptr(),Minimization::sort_per_energy); }
		best.set_target();
		unsigned int i(0);
		while(best.target_next() && i++<nsave){ best.get().save(&m_.system_param_); }
	} else {
		std::cerr<<"void VMCPSO::save(unsigned int const& nsave) : there is no data"<<std::endl;
	}
}
