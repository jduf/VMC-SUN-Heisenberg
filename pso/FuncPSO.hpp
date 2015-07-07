#ifndef DEF_FUNCPSO
#define DEF_FUNCPSO

#include "ParticleOnGrid.hpp"

class FuncPSO : public Swarm<ParticleOnGrid> {
	public:
		FuncPSO(unsigned int Nparticles, unsigned int maxiter, unsigned int Nfreedom, double cg, double cp);
		~FuncPSO() = default;
		void result();

		double f(Vector<double> const& x);

	protected:
		List<Measure> m_;

	private:
		bool evaluate(unsigned int const& p); 
};
#endif
