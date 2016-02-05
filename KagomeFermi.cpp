#include "KagomeFermi.hpp"

template<>
void KagomeFermi<double>::create(){
	compute_H();
	diagonalize(true);

	if(status_==2){
		for(unsigned int c(0);c<N_;c++){
			EVec_[c].set(n_,M_(c));
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = H_(i,j);
				}
			}
		}
	}
}

template<>
void KagomeFermi<std::complex<double> >::create(){
	compute_H();
	//select_eigenvectors();

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
}
