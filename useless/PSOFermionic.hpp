#ifndef DEF_PSOFERMIONIC
#define DEF_PSOFERMIONIC

#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"
#include "PSO.hpp"

class PSOFermionic : public PSO {
	public:
		PSOFermionic(Parseur& P,unsigned int Nfreedom,unsigned int N_MC);
		virtual ~PSOFermionic();

	private:
		CreateSystem CS_;
		unsigned int N_MC_;

		virtual double run(Vector<double> const& x);
};
#endif
