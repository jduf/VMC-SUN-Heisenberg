#include "VMCPSO.hpp"

VMCPSO::VMCPSO(Parseur& P, VMCMinimization const& vmcm):
	VMCMinimization(vmcm,"PSO"),
	Swarm<MCParticle>(P.get<unsigned int>("Nparticles"),P.get<unsigned int>("maxiter"),m_->Nfreedom_,P.get<double>("cg"),P.get<double>("cp"))
{
	for(unsigned int i(0);i<m_->Nfreedom_;i++){
		Particle::set_limit(i,0,m_->ps_[i].size());
	}
}

/*{Public methods*/
void VMCPSO::init(bool const& clear_particle_history, bool const& create_particle_history){
	set_time();

	std::cout<<"#######################"<<std::endl;
	std::cout<<"#new VMCPSO"<<std::endl;
	std::cout<<"#"<<get_filename()<<std::endl;
	std::cout<<"#contains "<<m_->samples_list_.size()<<" samples"<<std::endl;
	m_->pso_info_.title("New PSO run",'-');
	m_->pso_info_.item(get_filename());

	if(clear_particle_history){ 
		m_->pso_info_.item("clear history");
		for(unsigned int p(0);p<Nparticles_;p++){
			std::dynamic_pointer_cast<MCParticle>(particle_[p])->clear_history();
		}
	} else { m_->pso_info_.item("keep old history"); }

	std::shared_ptr<MCParticle> MCP;
	for(unsigned int i(0);i<Nparticles_;i++){
		MCP = std::dynamic_pointer_cast<MCParticle>(particle_[i]);
		MCP->set_ps(m_->ps_);
	}
	init_PSO(100); 
	if(clear_particle_history && create_particle_history && m_->samples_list_.size()){
		unsigned int size(0);
		while(m_->samples_list_.target_next()){
			if(m_->within_limit(m_->samples_list_.get().get_param())){ size++; }
		}
		unsigned int Npp(size/Nparticles_);
		unsigned int s;
		std::shared_ptr<MCParticle> MCP;
		for(unsigned int p(0);p<Nparticles_;p++){
			if(p == Nparticles_ - (size%Nparticles_) ){ Npp++; }
			MCP = std::dynamic_pointer_cast<MCParticle>(particle_[p]);
			s = 0;
			while( s!=Npp && m_->samples_list_.target_next()){
				if(m_->within_limit(m_->samples_list_.get().get_param())){
					MCP->add_to_history(m_->samples_list_.get_ptr());
					s++;
				}
			}
			MCP->select_new_best();
		}
		m_->pso_info_.item("set history");
		m_->pso_info_.item("each particle knows "+my::tostring(Npp)+" samples");
	} else { m_->pso_info_.item("start with empty history"); }
}

void VMCPSO::run(){
	std::string msg1("explore with "+my::tostring(Nparticles_)+" particles for "+my::tostring(maxiter_)+" steps");
	msg1 += " estimated time "+my::tostring(Nparticles_*(maxiter_+1)/omp_get_max_threads())+"s";
	std::cout<<"#"<<msg1<<std::flush;
	Time chrono;

	minimize();

	std::string msg2(" (done in "+my::tostring(chrono.elapsed())+"s)");
	std::cout<<msg2<<std::endl;
	m_->pso_info_.item(msg1+msg2);
}

void VMCPSO::plot() const {
	std::string filename(get_filename());
	IOFiles data(filename+".dat",true);
	m_->samples_list_.set_target();
	while(m_->samples_list_.target_next()){
		data<<m_->samples_list_.get().get_param()<<m_->samples_list_.get().get_S()->get_energy()<<IOFiles::endl;
	}
	m_->samples_list_.target_next();
	Gnuplot gp("./",filename);
	unsigned int N(m_->samples_list_.get().get_param().size());
	for(unsigned int i(0);i<N;i++){
		gp+=std::string(!i?"plot":"    ")+" '"+filename+".dat' u "+my::tostring(N+1)+":"+my::tostring(i+1)+":"+my::tostring(N+2)+" w xe notitle"+(i==N-1?"":",\\");
	}
	gp.save_file();
	gp.create_image(true);
}
/*}*/

/*{Private methods*/
bool VMCPSO::evaluate(unsigned int const& p){
	std::shared_ptr<MCParticle> MCP(std::dynamic_pointer_cast<MCParticle>(particle_[p]));
	std::shared_ptr<MCSim> sim(VMCMinimization::evaluate(MCP->get_param()));
	return (sim.get()?MCP->update(sim):false);
}
/*}*/
