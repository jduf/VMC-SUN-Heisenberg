#include "VMCPSO.hpp"

VMCPSO::VMCPSO(Parseur const& P, VMCMinimization const& vmcm, int const& which_symmetry):
	VMCMinimization(vmcm,"PSO"),
	Swarm<MCParticle>(P.get<unsigned int>("Nparticles"),P.get<unsigned int>("maxiter"),m_->dof_,P.get<double>("cg"),P.get<double>("cp"))
{
	for(unsigned int i(0);i<m_->dof_;i++){
		Particle::set_limit(i,0,m_->ps_[i].size());
	}
	if(which_symmetry){
		Vector<double> tmp(m_->dof_);
		std::vector<Matrix<int> > sym;
		CreateSystem cs(m_->s_);
		cs.init(&tmp,NULL);
		cs.get_wf_symmetries(sym);
		if(which_symmetry<0){
			for(unsigned int p(0);p<Nparticles_;p++){
				std::dynamic_pointer_cast<MCParticle>(particle_[p])->set_symmetry(sym[p%sym.size()]);
			}
		} else {
			for(unsigned int p(0);p<Nparticles_;p++){
				std::dynamic_pointer_cast<MCParticle>(particle_[p])->set_symmetry(sym[which_symmetry-1]);
			}
		}
		if(Nparticles_<sym.size()){
			std::cerr<<__PRETTY_FUNCTION__<<" : not enough particles with respect to the number of symmetries : "<<sym.size()<<std::endl;
		}
	}
}

/*{public methods*/
void VMCPSO::init_param_and_symmetry(Vector<double> const& param){
	if(m_->tmax_){
		set_time();

		std::cout<<"#######################"<<std::endl;
		std::string msg("VMCPSO minimize with param and symmetry");
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.title(msg,'-');

		std::cout<<"#"<<get_filename()<<std::endl;
		m_->info_.item(get_filename());

		msg="contains "+my::tostring(m_->samples_.size())+" samples";
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		Matrix<int> sym;
		std::shared_ptr<MCSim> sim(std::make_shared<MCSim>(param));
		List<MCSim>::Node* sample(NULL);
		if( m_->samples_.size() && !m_->samples_.find_in_sorted_list(sim,sample,MCSim::sort_for_merge)){
			/*!create a new sample*/
			msg = "param had never been tested";
			std::cout<<"#"<<msg<<std::endl;
			m_->info_.item(msg);

			evaluate_until_precision(param,0,1e-5,10);
		}
		if(m_->samples_.find_in_sorted_list(sim,sample,MCSim::sort_for_merge)){ sim = sample->get(); }
		else { std::cerr<<__PRETTY_FUNCTION__<<" bug"<<std::endl; }

		Time chrono;
		msg="initialize particles with given param and symmetry";
		m_->info_.item(msg);
		std::cout<<"#"<<msg<<std::endl;
		for(unsigned int p(0);p<Nparticles_;p++){
			std::shared_ptr<MCParticle> MCP(std::dynamic_pointer_cast<MCParticle>(particle_[p]));
			MCP->set_ps(m_->ps_);
			MCP->set_symmetry(sym);
		}
		init_PSO(100);
		for(unsigned int p(0);p<Nparticles_;p++){
			std::dynamic_pointer_cast<MCParticle>(particle_[p])->update(sim);
		}
		m_->effective_time_ = chrono.elapsed()*omp_get_max_threads()/Nparticles_;
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : tmax_ = 0"<<std::endl; }
}

void VMCPSO::init(bool const& clear_particle_history){
	if(m_->tmax_){
		set_time();

		std::cout<<"#######################"<<std::endl;
		std::string msg("VMCPSO");
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.title(msg,'-');

		std::cout<<"#"<<get_filename()<<std::endl;
		m_->info_.item(get_filename());

		msg="contains "+my::tostring(m_->samples_.size())+" samples";
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		if(m_->samples_.size()){
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
		for(unsigned int p(0);p<Nparticles_;p++){
			std::dynamic_pointer_cast<MCParticle>(particle_[p])->set_ps(m_->ps_);
		}
		init_PSO(100);
		m_->effective_time_ = chrono.elapsed()*omp_get_max_threads()/Nparticles_;
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : tmax_ = 0"<<std::endl; }
}

void VMCPSO::run(){
	if(m_->tmax_){
		std::string msg1("explore with "+my::tostring(Nparticles_)+" particles for "+my::tostring(maxiter_)+" steps,");
		msg1 += " estimated time "+my::tostring(1.1*Nparticles_*maxiter_*m_->effective_time_/omp_get_max_threads())+"s ";
		std::cout<<"#"<<msg1<<std::endl;
		Time chrono;

		minimize();

		std::string msg2("(done in "+my::tostring(chrono.elapsed())+"s)");
		std::cout<<msg2<<std::endl;
		m_->info_.item(msg1+msg2);
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : tmax_ = 0"<<std::endl; }
}
/*}*/

/*{private methods*/
bool VMCPSO::evaluate(unsigned int const& p){
	std::shared_ptr<MCParticle> MCP(std::dynamic_pointer_cast<MCParticle>(particle_[p]));
	std::shared_ptr<MCSim> sim(VMCMinimization::evaluate(MCP->get_param(),0));
	return (sim.get()?MCP->update(sim):false);
}
/*}*/
