#ifndef DEF_PSOFERMIONIC
#define DEF_PSOFERMIONIC

#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"
#include "PSO.hpp"

class MCSim{
	public:
		MCSim(Vector<double> const& param, Data<double> const& E): param_(param), E_(E) {}

		static unsigned int cmp_for_fuse(MCSim* list, MCSim* new_elem);
		static void fuse(MCSim* list, MCSim* new_elem);
		Data<double> const& get_energy(){ return E_; }

	private:
		Vector<double> param_;
		Data<double> E_;
};

class PSOFermionic : public PSO{
	public:
		PSOFermionic(Parseur& P);
		virtual ~PSOFermionic(){};

	private:
		virtual double f(Vector<double> const& x);
		void create();
		template<typename Type>
			double monte_carlo(CreateSystem& cs, Vector<double> const& x);

		unsigned int tmax_;
		Container parameters_;
		List<MCSim> all_results_;
};

template<typename Type>
double PSOFermionic::monte_carlo(CreateSystem& cs, Vector<double> const& x){
	MCSystem<Type>* S(NULL);
	if(cs.is_bosonic())
	{ S = new SystemBosonic<Type>(*dynamic_cast<const Bosonic<Type>*>(cs.get_system())); } 
	else 
	{ S = new SystemFermionic<Type>(*dynamic_cast<const Fermionic<Type>*>(cs.get_system())); }

	MonteCarlo<Type> sim(S,tmax_);
	sim.run();
	MCSim run_results(x,S->get_energy());

#pragma omp critical
	{
		all_results_.add_or_fuse_sort(&run_results, MCSim::cmp_for_fuse, MCSim::fuse);
	}

	delete S;
	return run_results.get_energy().get_dx();
}
#endif
