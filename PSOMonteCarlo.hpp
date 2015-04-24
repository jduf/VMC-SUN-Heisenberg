#ifndef DEF_PSOMONTECARLO
#define DEF_PSOMONTECARLO

#include "MonteCarlo.hpp"
#include "MCParticle.hpp"

class PSOMonteCarlo: public Swarm<MCParticle>{
	public:
		PSOMonteCarlo(Parseur* P);
		/*Default destructor*/
		virtual ~PSOMonteCarlo() = default;
		/*Forbids constructors*/
		PSOMonteCarlo() = delete;
		PSOMonteCarlo(PSOMonteCarlo const&) = delete;
		PSOMonteCarlo(PSOMonteCarlo&&) = delete;

		void complete_analysis(double tol);

		void plot() const;
		void print() const;
		void write(IOFiles& w) const;
		void load(std::string const& filename);

	private:
		bool is_better_x(unsigned int const& p);
		void create();

		unsigned int tmax_;
		Container system_;
		List<MCSim> all_results_;
};
#endif
