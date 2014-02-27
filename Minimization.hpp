#ifndef DEF_MINIMIZATION
#define DEF_MINIMIZATION

#include "Parseur.hpp"
#include "ParallelMonteCarlo.hpp"
#include "CreateSystem.hpp"

class Minimization{
	public:
		Minimization(Parseur& P);
		~Minimization();
	
		void min(double xmax);

	private:
		CreateSystem CS_;
		unsigned int tmax_;
		unsigned int nruns_;
		Write results_file_;

		double fx(double delta);
};
#endif
