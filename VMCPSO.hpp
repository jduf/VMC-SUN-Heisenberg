#ifndef DEF_VMCPSO
#define DEF_VMCPSO

#include "MCParticle.hpp"
#include "VMCMinimization.hpp"

class VMCPSO: public Swarm<MCParticle>, public VMCMinimization {
	public:
		VMCPSO(Parseur& P);
		/*!Default destructor*/
		virtual ~VMCPSO() = default;
		/*{Forbidden*/
		VMCPSO() = delete;
		VMCPSO(VMCPSO const&) = delete;
		VMCPSO(VMCPSO&&) = delete;
		VMCPSO& operator=(VMCPSO) = delete;
		/*}*/

		void init(bool const& clear_particle_history, bool const& create_particle_history);
		void refine(unsigned int const& Nrefine, double const& convergence_criterion, unsigned int const& tmax);
		void save_best(unsigned int const& nsave);
		void plot() const;
		void print() const;

	private:
		bool evaluate(unsigned int const& p);
		static bool sort_per_energy(MCSim const& a, MCSim const& b){ 
			return a.get_S()->get_energy().get_x()<b.get_S()->get_energy().get_x();
		};
};
#endif
