#include "MonteCarlo.hpp"

template<>
void MonteCarlo<double>::run(unsigned int const& thread){
	unsigned int i(1);
	double E_config(0);
	double ratio(0.0);
	Rand rnd(1e4,thread);
	do{
		S[thread].swap();
		ratio = S[thread].ratio(); 
		ratio *= ratio;
		if(ratio>1 || rnd.get() <ratio){
			S[thread].update();
			E_config = S[thread].compute_energy();
			E[thread] += E_config;
			sampling[thread].push_back(E_config);
			if(i < N_MC){ i++;}
			else {
				i = 1;
				test_convergence(thread);
			}
		}
	} while(keep_measuring);
#pragma omp critical
	{
		test_convergence(thread);
	}
}

template<>
void MonteCarlo<std::complex<double> >::run(unsigned int const& thread){
	unsigned int i(1);
	double E_config(0);
	double ratio(0.0);
	Rand rnd(1e4,thread);
	do{
		S[thread].swap();
		ratio = std::norm(S[thread].ratio());
		if(ratio>1 || rnd.get() <ratio){
			S[thread].update();
			E_config = S[thread].compute_energy();
			E[thread] += E_config;
			sampling[thread].push_back(E_config);
			if(i < N_MC){ i++;}
			else {
				i = 1;
				test_convergence(thread);
			}
		}
	} while(keep_measuring);
#pragma omp critical
	{
		test_convergence(thread);
	}
}
