#include "SquareFermi.hpp"

template<>
void SquareFermi<double>::create(){
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
void SquareFermi<std::complex<double> >::create(){
	std::cout<<"complex"<<std::endl;
	compute_H();
	diagonalize(false);
	select_eigenvectors();
	if(status_==1){
		for(unsigned int c(0);c<N_;c++){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = evec_(i,j);
				}
			}
		}
	}
	H_.set();
}
