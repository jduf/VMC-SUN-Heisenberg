#include "MonteCarlo.hpp"


template<>
void MonteCarlo<std::complex<double> >::compute_energy(){
	double E_step(0.0);
	for(unsigned int j(0);j<sts.row();j++){
		S->swap(sts(j,0),sts(j,1));
		E_step += std::real(S->ratio() * H(sts(j,0),sts(j,1)));
	}
	E += E_step;
	sampling.push_back(E_step);
}

template<>
void MonteCarlo<double>::compute_energy(){
	double E_step(0.0);
	for(unsigned int j(0);j<sts.row();j++){
		S->swap(sts(j,0),sts(j,1));
		E_step += S->ratio() * H(sts(j,0),sts(j,1));
	}
	E += E_step;
	sampling.push_back(E_step);
}

template<>
void MonteCarlo<double>::run(){
	srand(time(NULL)^(getpid()<<16));
	double ratio(0.0);
	unsigned int DCT(0);
	unsigned int i(0);
	do{
		S->swap();
		ratio = S->ratio(); 
		ratio *= ratio;
		if(ratio>1 || (double)rand()/RAND_MAX <ratio){
			S->update();
			if(DCT == 0){
				compute_energy();
				DCT = S->N_site;
				i++;
			}
			DCT--;
		}
	} while(binning_analyse(i));
	E /= (sampling.size() * S->N_site);
}

template<>
void MonteCarlo<std::complex<double> >::run(){
	srand(time(NULL)^(getpid()<<16));
	double ratio(0.0);
	unsigned int i(0);
	unsigned int DCT(0);
	do{
		S->swap();
		ratio = std::norm(S->ratio());
		if(ratio>1 || (double)rand()/RAND_MAX <ratio){
			S->update();
			if(DCT == 0){
				compute_energy();
				DCT = S->N_site;
				i++;
			}
			DCT--;
		}
	} while(binning_analyse(i));
	E /= (sampling.size() * S->N_site);
}
