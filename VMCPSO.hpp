#ifndef DEF_VMCPSO
#define DEF_VMCPSO

#include "MCParticle.hpp"
#include "VMCMinimization.hpp"

class VMCPSO: public VMCMinimization, public Swarm<MCParticle>{
	public:
		VMCPSO(Parseur& P, VMCMinimization const& vmcm);
		/*!Default destructor*/
		virtual ~VMCPSO() = default;
		/*{Forbidden*/
		VMCPSO() = delete;
		VMCPSO(VMCPSO const&) = delete;
		VMCPSO(VMCPSO&&) = delete;
		VMCPSO& operator=(VMCPSO) = delete;
		/*}*/

		void init(bool const& clear_particle_history, bool const& create_particle_history);
		void run();
		void plot() const;

	private:
		/*!Create param p then call VMCMinimization::evaluate(param)*/
		bool evaluate(unsigned int const& p);
};
#endif
