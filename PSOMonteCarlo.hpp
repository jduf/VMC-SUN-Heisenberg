#ifndef DEF_PSOMONTECARLO
#define DEF_PSOMONTECARLO

#include "MonteCarlo.hpp"
#include "MCParticle.hpp"

class PSOMonteCarlo: public Swarm<MCParticle>{
	public:
		PSOMonteCarlo(Parseur& P, IOFiles* in=NULL);
		/*!Default destructor*/
		virtual ~PSOMonteCarlo() = default;
		/*{Forbidden*/
		PSOMonteCarlo() = delete;
		PSOMonteCarlo(PSOMonteCarlo const&) = delete;
		PSOMonteCarlo(PSOMonteCarlo&&) = delete;
		PSOMonteCarlo& operator=(PSOMonteCarlo) = delete;
		/*}*/

		void refine(unsigned int const& Nrefine, double const& tol, unsigned int const& tmax);
		void complete_analysis(double const& tol);
		void plot() const;
		void print() const;

		void save();
		void create_particle_history(bool create);

	private:
		bool evaluate(unsigned int const& p);

		unsigned int tmax_;
		Container system_param_;
		List<MCSim> all_results_;
		IOFiles out_;

		std::string init_read(Parseur& P, IOFiles* in=NULL);
};
#endif
