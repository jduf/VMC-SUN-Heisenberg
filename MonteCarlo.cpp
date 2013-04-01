#include "MonteCarlo.hpp"

template<>
void MonteCarlo<double>::run(unsigned int const& thread){
	unsigned int i(0);
	double E_config(0);
	double ratio(0.0);
	do{
		S[thread].swap();
		ratio = S[thread].ratio(); 
		ratio *= ratio;
		if(ratio>1 || (double)rand()/RAND_MAX <ratio){
			S[thread].update();
			E_config = S[thread].compute_energy();
			E[thread] += E_config;
			sampling[thread].push_back(E_config);
			i++;
		}
	} while(binning_analysis(i,thread));
}

template<>
void MonteCarlo<std::complex<double> >::run(unsigned int const& thread){
	unsigned int i(0);
	double E_config(0);
	double ratio(0.0);
	do{
		S[thread].swap();
		ratio = std::norm(S[thread].ratio());
		if(ratio>1 || (double)rand()/RAND_MAX <ratio){
			S[thread].update();
			E_config = S[thread].compute_energy();
			E[thread] += E_config;
			sampling[thread].push_back(E_config);
			i++;
		}
	} while(binning_analysis(i,thread));
}
