#include"ChainFermi.hpp"

template<>
void ChainFermi<double>::create(){
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);
	lr_corr_.set(links_.row(),50,5,false);

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
void ChainFermi<std::complex<double> >::create(){
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);
	lr_corr_.set(links_.row(),50,5,false);

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
