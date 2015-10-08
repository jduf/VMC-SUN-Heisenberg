#include "KagomeDirac.hpp"

template<>
void KagomeDirac<double>::create(unsigned int const& which_observables){
	(void)(which_observables);
	compute_H();
	diagonalize(false);
	if(status_==2){
		std::cout<<"double"<<std::endl;
		for(unsigned int c(0);c<N_;c++){
			EVec_[c].set(n_,M_(c));
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = H_(i,j);
				}
			}
		}
		compute_H();
	}
}

template<>
void KagomeDirac<std::complex<double> >::create(unsigned int const& which_observables){
	(void)(which_observables);
	compute_H();
	select_eigenvectors();
	if(status_==2){
		for(unsigned int c(0);c<N_;c++){
			EVec_[c].set(n_,M_(c));
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = evec_(i,j);
				}
			}
		}
	}
	H_.set();
}
