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
		void save_best(unsigned int const& nsave);
		void plot() const;
		void print() const;
		void set_ps(unsigned int const& i, Vector<double> const& ps);

	private:
		bool within_limit(Vector<double> const& x);
		bool evaluate(unsigned int const& p);
};
#endif
