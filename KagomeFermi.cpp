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

	//compute_H();
	//BandStructure<std::complex<double> > bs(H_,Lx_,Ly_,spuc_,bc_);
	//Matrix<std::complex<double> > eval(n_,3);
	//Vector<double> tmp_eval;
	//Lapack<std::complex<double> >(H_,true,'H').eigensystem(tmp_eval,false);
	//for(unsigned int i(0);i<n_;i++){
		//eval(i,0) = tmp_eval(i);
	//}
	//std::cout<<eval<<std::endl;
	//bs.diagonalize_everything(evec_,eval);

}
