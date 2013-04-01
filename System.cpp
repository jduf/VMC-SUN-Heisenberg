#include "System.hpp"

template<>
double System<double>::compute_energy(){
	double E_step(0.0);
	for(unsigned int j(0);j<sts.row();j++){
		swap(sts(j,0),sts(j,1));
		E_step += ratio() * H(sts(j,0),sts(j,1));
	}
	return E_step;
}

template<>
double System<std::complex<double> >::compute_energy(){
	double E_step(0.0);
	for(unsigned int j(0);j<sts.row();j++){
		swap(sts(j,0),sts(j,1));
		E_step += std::real(ratio() * H(sts(j,0),sts(j,1)));
	}
	return E_step;
}
