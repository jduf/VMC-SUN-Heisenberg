#ifndef DEF_PSOFERMIONIC
#define DEF_PSOFERMIONIC

#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"
#include "PSO.hpp"
#include "List.hpp"

class MCSim {
	public:
		MCSim(Vector<double> const& param): param_(param) { 
//#pragma omp critical(all_results_)
			//{
				//std::cout<<"create "<<param_<<" "<<this<<" "<<S_.get()<<std::endl;
			//}
		}
		~MCSim(){
//#pragma omp critical(all_results_)
			//{
				//std::cout<<"destroy "<<this<<std::endl; 
			//}
		}

		static unsigned int cmp_for_fuse(MCSim const& list_elem, MCSim const& new_elem);
		static void fuse(MCSim& list_elem, MCSim& new_elem);

		Vector<double> const& get_param() const { return param_; }
		std::unique_ptr<MCSystem> const& get_S() const { return S_; }
		void create_S(Container* C);
		void copy_S(std::unique_ptr<MCSystem> const& S);

		void print(std::ostream& flux) const { flux<<param_<<S_->get_energy(); }

	private:
		Vector<double> param_;
		std::unique_ptr<MCSystem> S_;
};

class MCParticle: public Particle{
	public:
		MCParticle(){}
		~MCParticle(){}

		void move(Vector<double> const& bx_all);
		bool update(std::shared_ptr<MCSim> const& new_elem);

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

		unsigned int tmax_;
		Container system_;
		List<MCSim> all_results_;
};

std::ostream& operator<<(std::ostream& flux, MCSim const& mcsim);
std::ostream& operator<<(std::ostream& flux, PSOFermionic const& pso);
#endif
