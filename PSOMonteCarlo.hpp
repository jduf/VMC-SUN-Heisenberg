#ifndef DEF_PSOMONTECARLO
#define DEF_PSOMONTECARLO

#include "MonteCarlo.hpp"
#include "MCParticle.hpp"

class PSOMonteCarlo: public Swarm<MCParticle>{
	public:
		PSOMonteCarlo(Parseur* P);
		/*!Default destructor*/
		virtual ~PSOMonteCarlo() = default;
		/*{Forbids constructors*/
		PSOMonteCarlo() = delete;
		PSOMonteCarlo(PSOMonteCarlo const&) = delete;
		PSOMonteCarlo(PSOMonteCarlo&&) = delete;
		PSOMonteCarlo& operator=(PSOMonteCarlo) = delete;
		/*}*/

		void refine(unsigned int const& Nrefine, double const& tol, unsigned int const& tmax);
		void complete_analysis(double const& tol);
		void plot() const;
		void print() const;

		void write(IOFiles& w) const;
		void read(IOFiles& w, bool create_particle_history);

	private:
		bool evaluate(unsigned int const& p);

		unsigned int tmax_;
		Container system_param_;
		List<MCSim> all_results_;
};
#endif
