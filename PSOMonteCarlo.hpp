#ifndef DEF_PSOMONTECARLO
#define DEF_PSOMONTECARLO

#include "MonteCarlo.hpp"
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
		void complete_analysis(double const& tol);
		void refine(unsigned int const& Nrefine, double const& tol, unsigned int const& tmax);
		void save() const;
		void save(unsigned int const& nsave);
		void plot() const;
		void print() const;

	private:
		bool evaluate(unsigned int const& p);

		unsigned int tmax_;
		Container system_param_;
		List<MCSim> all_results_;
		RST pso_info_;

		std::string get_filename() const;
};
#endif
