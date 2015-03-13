#ifndef DEF_PSOFERMIONIC
#define DEF_PSOFERMIONIC

#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"
#include "PSO.hpp"
#include "List.hpp"

class MCSim{
	public:
		MCSim(Vector<double> const& param, Data<double> const& E): param_(param), E_(E), N_(1) {}

		static unsigned int cmp_for_fuse(MCSim* list, MCSim* new_elem);
		static void fuse(MCSim* list, MCSim* new_elem);
		Vector<double> const& get_param(){ return param_; }
		Data<double> const& get_energy(){ return E_; }
		unsigned int const& get_N(){ return N_; }

		void print(std::ostream& flux) const { flux<<param_<< E_<<" ("<<N_<<")"<<std::endl;; }

	private:
		Vector<double> param_;
		Data<double> E_;
		unsigned int N_;
};

std::ostream& operator<<(std::ostream& flux, MCSim const& mcsim);

class PSOFermionic : public PSO{
	public:
		PSOFermionic(Parseur* P);
		virtual ~PSOFermionic(){}

		void print(){ std::cout<<all_results_<<std::endl; }
		void plot(){
			IOFiles data("data.dat",true);
			for(unsigned int i(0);i<all_results_.size();i++){
				data<<all_results_[i].get_param()<<" "<<all_results_[i].get_energy()<<IOFiles::endl;
			}
			Gnuplot gp("./","test");
			gp+="plot 'data.dat' u 3:1,\\";
			gp+="     'data.dat' u 3:2";
			gp.save_file();
			gp.create_image(true);
		}

	private:
		virtual double f(Vector<double> const& x);
		void create();
		template<typename Type>
			double monte_carlo(CreateSystem& cs, Vector<double> const& x);

		unsigned int tmax_;
		Container system_;
		List<MCSim> all_results_;
};

template<typename Type>
double PSOFermionic::monte_carlo(CreateSystem& cs, Vector<double> const& x){
	MCSystem<Type>* S(NULL);
	if(cs.is_bosonic())
	{ S = new SystemBosonic<Type>
		(*dynamic_cast<const Bosonic<Type>*>(cs.get_system())); } 
	else 
	{ S = new SystemFermionic<Type>
		(*dynamic_cast<const Fermionic<Type>*>(cs.get_system())); }

	MonteCarlo<Type> sim(S,tmax_);
	if(S->get_status() == 0){
		sim.thermalize(1e6);
		sim.run();
		sim.complete_analysis(1e-5);
		MCSim* run_results(new MCSim(x,S->get_energy()));

#pragma omp critical(add_new_result_to_list)
		{
			all_results_.add_or_fuse_sort(run_results, MCSim::cmp_for_fuse, MCSim::fuse);
		}

		delete S;
		return run_results->get_energy().get_x();
	} else {
		std::cerr<<"double PSOFermionic::monte_carlo(CreateSystem& cs, Vector<double> const& x) : no value for x="<<x<<std::endl;
		delete S;
		return 0.0;
	}
}
#endif
