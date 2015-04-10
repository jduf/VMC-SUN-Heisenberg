#ifndef DEF_PSOFERMIONIC
#define DEF_PSOFERMIONIC

#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"
#include "PSO.hpp"
#include "List.hpp"

class MCSim {
	public:
		MCSim(Vector<double> const& param, Data<double> const& E): param_(param), E_(E), N_(1) {}

		static unsigned int cmp_for_fuse(MCSim const& list_elem, MCSim const& new_elem);
		static void fuse(MCSim& list_elem, MCSim& new_elem);

		Vector<double> const& get_param() const { return param_; }
		Data<double> const& get_energy() const { return E_; }
		unsigned int const& get_N() const { return N_; }

		void print(std::ostream& flux) const { flux<<param_<< E_<<" ("<<N_<<")"<<std::endl;; }

	private:
		Vector<double> param_;
		Data<double> E_;
		unsigned int N_;
};

class MCParticle: public Particle{
	public:
		MCParticle(){}
		~MCParticle(){}

		void move(Vector<double> const& bx_all);
		void update_particle_history(std::shared_ptr<MCSim>& new_elem);

		static unsigned int pos_iter;
		static Vector<double> pos;

		void print(std::ostream& flux) const;

	private:
		List<MCSim> history_;
		static void fuse(MCSim& list_elem, MCSim& new_elem);
};

class PSOFermionic: public Swarm<MCParticle>{
	public:
		PSOFermionic(Parseur* P);
		virtual ~PSOFermionic(){}

		void plot();
		void print(std::ostream& flux) const;

	private:
		bool is_better_x(unsigned int const& p);
		void create();
		template<typename Type>
			bool monte_carlo(CreateSystem& cs, unsigned int const& p);

		unsigned int tmax_;
		Container system_;
		List<MCSim> all_results_;
};

template<typename Type>
bool PSOFermionic::monte_carlo(CreateSystem& cs, unsigned int const& p){
	std::shared_ptr<MCParticle> P(std::dynamic_pointer_cast<MCParticle>(p_[p]));

	MCSystem<Type>* S;
	if(cs.is_bosonic()) {S = new SystemBosonic<Type>
		(*dynamic_cast<const Bosonic<Type>*>(cs.get_system()));}
	else { S = new SystemFermionic<Type>
		(*dynamic_cast<const Fermionic<Type>*>(cs.get_system()));}
	std::shared_ptr<MCSim> rr;
	double local_e(0.0);
	bool local_improvement(false);

	MonteCarlo<Type> sim(S,tmax_);
	if(S->get_status() == 0){
		sim.thermalize(1e6);
		do {
			sim.run();
			sim.complete_analysis(1e-5);
			rr = std::make_shared<MCSim>(P->get_x(),S->get_energy());
			local_e = rr->get_energy().get_x();
#pragma omp critical(add_new_result_to_list)
			{
				std::cout<<"x "<<P->get_x()<<std::endl;
				std::cout<<"e "<<rr->get_energy()<<std::endl;
				if(all_results_.find_sorted(rr,MCSim::cmp_for_fuse)){ all_results_.fuse_with_move(rr, MCSim::fuse); }
				else { all_results_.add_after_free(rr); }
				std::cout<<"e "<<rr->get_energy()<<std::endl;
			}
			P->update_particle_history(rr);
			local_improvement = rr->get_energy().get_x()<local_e;
		} while (local_improvement);
	} else {
		std::cerr<<"double PSOFermionic::monte_carlo(CreateSystem& cs, Vector<double> const& x) : no value for x="<<p_[p]->get_x()<<std::endl;
		return true;
	}

	delete S;
	P->update_particle_history(rr);
	if(local_e < P->get_fbx()){ 
		std::cout<<"update best local position"<<std::endl;
		P->set_best(P->get_x(),local_e);
		return true;
	} else {
		return false;
	}
}

std::ostream& operator<<(std::ostream& flux, MCSim const& mcsim);
std::ostream& operator<<(std::ostream& flux, PSOFermionic const& pso);
#endif
