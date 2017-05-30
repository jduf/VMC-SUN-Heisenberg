#ifndef DEF_VMCPSO
#define DEF_VMCPSO

#include "MCParticle.hpp"
#include "VMCMinimization.hpp"

class VMCPSO: public VMCMinimization, public Swarm<MCParticle>{
	public:
		VMCPSO(Parseur const& P, VMCMinimization const& vmcm);
		/*!Default destructor*/
		virtual ~VMCPSO() = default;
		/*{Forbidden*/
		VMCPSO() = delete;
		VMCPSO(VMCPSO const&) = delete;
		VMCPSO(VMCPSO&&) = delete;
		VMCPSO& operator=(VMCPSO) = delete;
		/*}*/

		void run(bool const& clear_particle_history, double const& dEoE, unsigned int const& tmax, unsigned int const& maxiter, std::string const& save_in);

	private:
		unsigned int tmax_;
		/*!Create param p then call VMCMinimization::evaluate(param)*/
		bool evaluate(unsigned int const& p);
};
#endif
