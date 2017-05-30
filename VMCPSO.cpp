#include "VMCPSO.hpp"

VMCPSO::VMCPSO(Parseur const& P, VMCMinimization const& vmcm):
	VMCMinimization(vmcm,"PSO"),
	Swarm<MCParticle>(P.get<unsigned int>("Nparticles"),P.get<unsigned int>("maxsteps"),m_->dof_,P.check_get("cg",2.1),P.check_get("cp",2.1))
{
	for(unsigned int i(0);i<m_->dof_;i++){
		Particle::set_limit(i,0,m_->ps_[i].size());
	}
}

/*{public methods*/
void VMCPSO::run(bool const& clear_particle_history, double const& dEoE, unsigned int const& tmax, unsigned int const& maxiter, std::string const& save_in){
	tmax_ = tmax;
	if(tmax_){
		set_time();
		std::cout<<"#######################"<<std::endl;
		std::string msg("VMCPSO");
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.title(msg,'-');

		std::cout<<"#"<<get_filename()<<std::endl;
		m_->info_.item(get_filename());
		msg = "contains "+my::tostring(m_->samples_.size())+" samples";
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		if(m_->samples_.size()){
			if(clear_particle_history){
				msg = "clear particles' history";
				m_->info_.item(msg);
				std::cout<<"#"<<msg<<std::endl;
				for(unsigned int p(0);p<Nparticles_;p++){
					std::dynamic_pointer_cast<MCParticle>(particle_[p])->clear_history();
				}
			} else {
				msg = "keep old particles' history";
				m_->info_.item(msg);
				std::cout<<"#"<<msg<<std::endl;
			}
		}

		Time chrono;
		msg = "initializing each particles during "+my::tostring(tmax_)+"s";
		m_->info_.item(msg);
		std::cout<<"#"<<msg<<std::endl;
		for(unsigned int p(0);p<Nparticles_;p++){
			std::dynamic_pointer_cast<MCParticle>(particle_[p])->set_ps(m_->ps_);
		}
		init_PSO(100);
		m_->effective_time_ = chrono.elapsed()*omp_get_max_threads()/Nparticles_;

		std::string msg1("exploring with "+my::tostring(Nparticles_)+" particles for "+my::tostring(maxsteps_)+" steps,");
		msg1 += " estimated time "+my::tostring(1.1*Nparticles_*maxsteps_*m_->effective_time_/omp_get_max_threads())+"s ";
		std::cout<<"#"<<msg1<<std::endl;

		chrono.set();
		minimize();

		std::string msg2("(done in "+my::tostring(chrono.elapsed())+"s)");
		std::cout<<msg2<<std::endl;
		m_->info_.item(msg1+msg2);

		improve_bad_samples(dEoE,tmax_,save_in); 
		refine(0,dEoE,maxiter,10*tmax_,30,save_in);
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : tmax_ = 0"<<std::endl; }
}
/*}*/

/*{private methods*/
bool VMCPSO::evaluate(unsigned int const& p){
	std::shared_ptr<MCParticle> MCP(std::dynamic_pointer_cast<MCParticle>(particle_[p]));
	std::shared_ptr<MCSim> sim(VMCMinimization::evaluate(MCP->get_param(),0,tmax_));
	return (sim.get()?MCP->update(sim):false);
}
/*}*/
