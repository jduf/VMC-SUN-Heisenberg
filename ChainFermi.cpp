#include"ChainFermi.hpp"

template<>
void ChainFermi<double>::create(unsigned int const& which_observables){
	(void)(which_observables);
	compute_H();
	diagonalize(true);
	if(status_==1){
		for(unsigned int c(0);c<N_;c++){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = H_(i,j);
				}
			}
		}
	}
}

template<>
void ChainFermi<std::complex<double> >::create(unsigned int const& which_observables){
	(void)(which_observables);
	compute_H();
	diagonalize(false);
	for(unsigned int c(0);c<N_;c++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_(c);j++){
				/*!Use the eigenvector k1*/
				EVec_[c](i,j) = evec_(i,j);
			}
		}
	}
}
