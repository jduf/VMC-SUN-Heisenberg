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
		Write file_;
		Container param_;
		Container param_text_;
		unsigned int nthreads;
		unsigned int t_max;

		double fx(double delta);
		void save(double delta, Container const& result);
};
#endif
