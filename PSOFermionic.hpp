#ifndef DEF_PSOFERMIONIC
#define DEF_PSOFERMIONIC

#include "MonteCarlo.hpp"
#include "ChainTrimerized.hpp"
#include "PSO.hpp"

class PSOFermionic : public PSO {
	public:
		PSOFermionic(Parseur& P,unsigned int Nfreedom);
		virtual ~PSOFermionic();

	private:
		virtual double run(Vector<double> const& x);

		Write results_;
		ChainTrimerized CS_;
};
#endif
