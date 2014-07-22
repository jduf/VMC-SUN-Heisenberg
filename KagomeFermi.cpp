#include "KagomeFermi.hpp"

template<>
void KagomeFermi<double>::create(){
	std::cout<<"double"<<std::endl;
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);

	compute_T();
	diagonalize_T();

	if(!degenerate_){
		for(unsigned int c(0);c<N_;c++){
			EVec_[c].set(n_,M_(c));
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = T_(i,j);
				}
			}
		}
	}
}

template<>
void KagomeFermi<std::complex<double> >::create(){
	std::cout<<"complex"<<std::endl;
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);

	compute_T();
	H_ = T_;
	compute_TxTy();
	compute_band_structure();
	plot_band_structure();
	select();


	//for(unsigned int c(0);c<N_;c++){
		//EVec_[c].set(n_,M_(c));
		//for(unsigned int i(0);i<n_;i++){
			//for(unsigned int j(0);j<M_(c)-1;j++){
				//EVec_[c](i,j) = evec_(i,j);
			//}
			//EVec_[c](i,M_(c)-1) = evec_(i,M_(c)-1+c);
		//}
		//double Px(px_(M_(c)-1));
		//double Py(py_(M_(c)-1));
		//double e(  e_(M_(c)-1));
		//for(unsigned int j(0);j<M_(c)-1;j++){
			//Px += (are_equal(std::abs(px_(j)),M_PI,1e-8)?0:px_(j));
			//Py += (are_equal(std::abs(py_(j)),M_PI,1e-8)?0:px_(j));
			//e += e_(j);
		//}
		//std::cout<<e<<" "<<Px<<" "<<Py<<std::endl;
	//}
	T_.set();
}
