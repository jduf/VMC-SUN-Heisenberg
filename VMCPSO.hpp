#ifndef DEF_VMCPSO
#define DEF_VMCPSO

#include "MCParticle.hpp"
#include "VMCMinimization.hpp"

class VMCPSO: public VMCMinimization, public Swarm<MCParticle>{
	public:
		VMCPSO(Parseur const& P, VMCMinimization const& vmcm, int const& set_symmetry);
		/*!Default destructor*/
		virtual ~VMCPSO() = default;
		/*{Forbidden*/
		VMCPSO() = delete;
		VMCPSO(VMCPSO const&) = delete;
		VMCPSO(VMCPSO&&) = delete;
		VMCPSO& operator=(VMCPSO) = delete;
		/*}*/

		void init_param_and_symmetry(Vector<double> const& param);
		void init(bool const& clear_particle_history);
		void run();

	private:
		/*!Create param p then call VMCMinimization::evaluate(param)*/
		bool evaluate(unsigned int const& p);
};
#endif
