#ifndef DEF_MINIMIZATION
#define DEF_MINIMIZATION

#include "Parseur.hpp"
#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"

class Minimization{
	public:
		Minimization(Parseur& P);
		~Minimization();
	
		void min(double xmax);

	private:
		unsigned int nthreads;
		unsigned int t_max;
		CreateSystem CS_;

		double fx(double delta);
};
#endif
