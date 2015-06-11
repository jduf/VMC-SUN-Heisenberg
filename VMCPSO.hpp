#ifndef DEF_VMCPSO
#define DEF_VMCPSO

#include "MonteCarlo.hpp"
#include "MCParticle.hpp"

class VMCPSO: public Swarm<MCParticle>{
	public:
		VMCPSO(Parseur* P);
		/*!Default destructor*/
		virtual ~VMCPSO() = default;
		/*{Forbids constructors*/
		VMCPSO() = delete;
		VMCPSO(VMCPSO const&) = delete;
		VMCPSO(VMCPSO&&) = delete;
		VMCPSO& operator=(VMCPSO) = delete;
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
