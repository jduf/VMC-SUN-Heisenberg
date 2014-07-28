#include "KagomeFermi.hpp"

template<>
void KagomeFermi<double>::create(){
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);

	compute_H();
	diagonalize_H(H_);

	if(!degenerate_){
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
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);

	compute_H();
	compute_TxTy();
	compute_band_structure();
	//plot_band_structure();
	select_eigenvectors();

	if(!degenerate_){
		for(unsigned int c(0);c<N_;c++){
			EVec_[c].set(n_,M_(c));
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = evec_(i,select_[c](j));
				}
			}
		}
	}
	H_.set();
}
