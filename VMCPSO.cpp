#include "VMCPSO.hpp"

VMCPSO::VMCPSO(Parseur& P, VMCMinimization const& vmcm):
	VMCMinimization(vmcm,"PSO"),
	Swarm<MCParticle>(P.get<unsigned int>("Nparticles"),P.get<unsigned int>("maxiter"),m_->Nfreedom_,P.get<double>("cg"),P.get<double>("cp"))
{
	for(unsigned int i(0);i<m_->Nfreedom_;i++){
		Particle::set_limit(i,0,m_->ps_[i].size());
	}
}

/*{public methods*/
void VMCPSO::init(bool const& clear_particle_history, bool const& create_particle_history){
	set_time();

	std::cout<<"#######################"<<std::endl;
	std::string msg("new VMCPSO");
	std::cout<<"#"<<msg<<std::endl;
	m_->pso_info_.title(msg,'-');

	std::cout<<"#"<<get_filename()<<std::endl;
	m_->pso_info_.item(get_filename());

	msg="contains "+my::tostring(m_->samples_list_.size())+" samples";
	std::cout<<"#"<<msg<<std::endl;
	m_->pso_info_.item(msg);

	if(clear_particle_history){ 
		msg="clear history";
		m_->pso_info_.item(msg); 
		std::cout<<"#"<<msg<<std::endl;
		for(unsigned int p(0);p<Nparticles_;p++){
			std::dynamic_pointer_cast<MCParticle>(particle_[p])->clear_history();
		}
	} else {
		msg="keep old history";
		m_->pso_info_.item(msg); 
		std::cout<<"#"<<msg<<std::endl;
	}

	msg="initialize particles";
	m_->pso_info_.item(msg); 
	std::cout<<"#"<<msg<<std::endl;
#pragma omp parallel for
	for(unsigned int i=0;i<Nparticles_;i++){
		std::shared_ptr<MCParticle> MCP;
		MCP = std::dynamic_pointer_cast<MCParticle>(particle_[i]);
		MCP->set_ps(m_->ps_);
	}
	
	Time chrono;
	init_PSO(100); 
	m_->effective_time_ = chrono.elapsed()*omp_get_max_threads()/Nparticles_;

	if(clear_particle_history && create_particle_history && m_->samples_list_.size()){
		msg="create particle history";
		std::cout<<"#"<<msg<<std::flush;

		unsigned int size(0);
		m_->samples_list_.set_target();
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
		std::string msg2(" (each particle knows "+my::tostring(Npp)+" samples)");
		m_->pso_info_.item(msg+msg2);
		std::cout<<msg2;
	} else {
		msg = "start with empty history";
		m_->pso_info_.item(msg); 
		std::cout<<"#"<<msg<<std::endl;
	}
}

void VMCPSO::run(){
	std::string msg1("explore with "+my::tostring(Nparticles_)+" particles for "+my::tostring(maxiter_)+" steps");
	msg1 += " estimated time "+my::tostring(1.1*Nparticles_*maxiter_*m_->effective_time_/omp_get_max_threads())+"s";
	std::cout<<"#"<<msg1<<std::flush;
	Time chrono;

	minimize();

	std::string msg2(" (done in "+my::tostring(chrono.elapsed())+"s)");
	std::cout<<msg2<<std::endl;
	m_->pso_info_.item(msg1+msg2);
}
/*}*/

/*{private methods*/
bool VMCPSO::evaluate(unsigned int const& p){
	std::shared_ptr<MCParticle> particle(std::dynamic_pointer_cast<MCParticle>(particle_[p]));
	std::shared_ptr<MCSim> sim(VMCMinimization::evaluate(particle->get_param()));
	return (sim.get()?particle->update(sim):false);
}
/*}*/
