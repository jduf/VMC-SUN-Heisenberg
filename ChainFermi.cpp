#include"ChainFermi.hpp"

template<>
void ChainFermi<double>::create(){
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);
	lr_corr_.set(links_.row(),50,5,false);

	compute_H();
	diagonalize_H(H_);
	if(!degenerate_){
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
	std::cout<<"complex"<<std::endl;
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);
	lr_corr_.set(links_.row(),50,5,false);

	///*{Description*/
	//compute_H();
	//select_eigenvectors(M_(0));
	//if(!degenerate_){
		//for(unsigned int c(0);c<N_;c++){
			//for(unsigned int i(0);i<n_;i++){
				//for(unsigned int j(0);j<M_(c);j++){
					//EVec_[c](i,j) = evec_(i,j);
				//}
			//}
			//std::cout<<(EVec_[c].trans_conj()*EVec_[c]).chop()<<std::endl;
		//}
	//}
	//H_.set();
	///*}*/

	/*{Description*/
	compute_H();
	diagonalize_H(H_);
	for(unsigned int c(0);c<N_;c++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_(c)-1;j++){
				EVec_[c](i,j) = H_(i,j);
			}
		}
	}
	compute_H();
	select_eigenvectors(M_(0));
	for(unsigned int c(0);c<N_;c++){
		for(unsigned int i(0);i<n_;i++){
			EVec_[c](i,M_(c)-1) = evec_(i,M_(c)-1);
		}
		//std::cout<<(EVec_[c].trans_conj()*EVec_[c]).chop()<<std::endl;
	}
	degenerate_ = false;
	/*}*/
	std::cerr<<"it is possible to use only double instead of complex double"<<std::endl;
}
