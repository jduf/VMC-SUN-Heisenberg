#ifndef DEF_PSOMONTECARLO
#define DEF_PSOMONTECARLO

#include "MCParticle.hpp"

class PSOMonteCarlo: public Swarm<MCParticle>{
	public:
		PSOMonteCarlo(Parseur& P);
		/*!Default destructor*/
		virtual ~PSOMonteCarlo() = default;
		/*{Forbidden*/
		PSOMonteCarlo() = delete;
		PSOMonteCarlo(PSOMonteCarlo const&) = delete;
		PSOMonteCarlo(PSOMonteCarlo&&) = delete;
		PSOMonteCarlo& operator=(PSOMonteCarlo) = delete;
		/*}*/

		void init(bool const& clear_particle_history, bool const& create_particle_history);
		void complete_analysis(double const& converged_criterion);
		void refine(unsigned int const& Nrefine, double const& converged_criterion, unsigned int const& tmax);
		void save() const;
		void save(unsigned int const& nsave);
		void plot() const;
		void print() const;

	private:
		bool evaluate(unsigned int const& p);

		unsigned int tmax_;
		Container system_param_;
		List<MCSim> all_results_;
		std::string basename_;
		std::string time_;
		RST pso_info_;

		std::string get_filename() const { return basename_+"_"+time_; }

		static bool sort_per_energy(MCSim const& a, MCSim const& b){ 
			return a.get_S()->get_energy().get_x()<b.get_S()->get_energy().get_x();
		};
};
#endif
