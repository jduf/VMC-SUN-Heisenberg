#include "VMCPSO.hpp"

VMCPSO::VMCPSO(Parseur const& P, VMCMinimization const& vmcm, bool set_symmetry):
	VMCMinimization(vmcm,"PSO"),
	Swarm<MCParticle>(P.get<unsigned int>("Nparticles"),P.get<unsigned int>("maxiter"),m_->dof_,P.get<double>("cg"),P.get<double>("cp"))
{
	for(unsigned int i(0);i<m_->dof_;i++){
		Particle::set_limit(i,0,m_->ps_[i].size());
	}
	if(set_symmetry){
		Vector<double> tmp(m_->dof_);
		std::vector<Matrix<int> > sym;
		CreateSystem cs(m_->s_);
		cs.init(&tmp,NULL);
		cs.get_wf_symmetries(sym);
		for(unsigned int i(0);i<Nparticles_;i++){
			std::shared_ptr<MCParticle> MCP;
			MCP = std::dynamic_pointer_cast<MCParticle>(particle_[i]);
			MCP->set_symmetry(sym[i%sym.size()]);
		}
		if(Nparticles_<sym.size()){
			std::cerr<<__PRETTY_FUNCTION__<<" : not enough particles with respect to the number of symmetries : "<<sym.size()<<std::endl;
		}
	}
}

/*{public methods*/
void VMCPSO::init(bool const& clear_particle_history){
	if(m_->tmax_){
		set_time();

		std::cout<<"#######################"<<std::endl;
		std::string msg("VMCPSO");
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.title(msg,'-');

		std::cout<<"#"<<get_filename()<<std::endl;
		m_->info_.item(get_filename());

		msg="contains "+my::tostring(m_->samples_list_.size())+" samples";
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		if(m_->samples_list_.size()){
			if(clear_particle_history){
				msg="clear particles' history";
				m_->info_.item(msg);
				std::cout<<"#"<<msg<<std::endl;
				for(unsigned int p(0);p<Nparticles_;p++){
					std::dynamic_pointer_cast<MCParticle>(particle_[p])->clear_history();
				}
			} else {
				msg="keep old particles' history";
				m_->info_.item(msg);
				std::cout<<"#"<<msg<<std::endl;
			}
		}

		Time chrono;
		msg="initialize particles";
		m_->info_.item(msg);
		std::cout<<"#"<<msg<<std::endl;
#pragma omp parallel for
		for(unsigned int i=0;i<Nparticles_;i++){
			std::dynamic_pointer_cast<MCParticle>(particle_[i])->set_ps(m_->ps_);
		}
		init_PSO(100);
		m_->effective_time_ = chrono.elapsed()*omp_get_max_threads()/Nparticles_;
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : tmax_ = 0"<<std::endl; }
}

void VMCPSO::run(){
	if(m_->tmax_){
		std::string msg1("explore with "+my::tostring(Nparticles_)+" particles for "+my::tostring(maxiter_)+" steps,");
		msg1 += " estimated time "+my::tostring(1.1*Nparticles_*maxiter_*m_->effective_time_/omp_get_max_threads())+"s";
		std::cout<<"#"<<msg1<<std::flush;
		Time chrono;

		minimize();

		std::string msg2(" (done in "+my::tostring(chrono.elapsed())+"s)");
		std::cout<<msg2<<std::endl;
		m_->info_.item(msg1+msg2);
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : tmax_ = 0"<<std::endl; }
}
/*}*/

/*{private methods*/
bool VMCPSO::evaluate(unsigned int const& p){
	std::shared_ptr<MCParticle> particle(std::dynamic_pointer_cast<MCParticle>(particle_[p]));
	std::shared_ptr<MCSim> sim(VMCMinimization::evaluate(particle->get_param(),0));
	return (sim.get()?particle->update(sim):false);
}
/*}*/
