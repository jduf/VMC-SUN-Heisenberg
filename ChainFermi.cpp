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
//
		///*suprisingly not too bad*/
		//Vector<double> v(n_);
		//for(unsigned int c(0);c<N_;c++){
			//switch(c){
				//case 0:
					//{
						//for(unsigned int i(0);i<n_;i++){
							//v(i) = H_(i,M_(c)) + H_(i,M_(c)-1);
						//}
					//}break;
				//case 1:
					//{
						//for(unsigned int i(0);i<n_;i++){
							//v(i) = H_(i,M_(c));
						//}
					//}break;
				//case 2:
					//{
						//for(unsigned int i(0);i<n_;i++){
							//v(i) = H_(i,M_(c)-1);
						//}
					//}break;
			//}
			//double norm(v.norm());
			//for(unsigned int i(0);i<n_;i++){ EVec_[c](i,M_(c)-1) = v(i)/norm; }
			//std::cout<<(EVec_[c].trans_conj()*EVec_[c]).chop()<<std::endl;
			//std::cout<<"..................."<<std::endl;
			//std::cout<<EVec_[c].chop()<<std::endl;
			//std::cout<<"..................."<<std::endl;
		//}
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
	//H_.set();
	//if(!degenerate_){
	//for(unsigned int c(0);c<N_;c++){
	//for(unsigned int i(0);i<n_;i++){
	//for(unsigned int j(0);j<M_(c);j++){
	//EVec_[c](i,j) = evec_(i,j);
	//}
	//}
	//}
	//for(unsigned int c(0);c<N_;c++){
	//Vector<std::complex<double> > v(n_);
	//for(unsigned int i(0);i<n_;i++){
	//v(i) = evec_(i,M_(c)-1);
	//v(i)+= evec_(i,M_(c));
	//}
	//double norm(v.norm());
	//for(unsigned int i(0);i<n_;i++){ EVec_[c](i,M_(c)-1) = v(i)/norm; }
	//std::cout<<(EVec_[c].trans_conj()*EVec_[c]).chop()<<std::endl;
	//}
	//}
	///*}*/

	///*{Description*/
	//compute_H();
	//diagonalize_H(H_);
	//degenerate_=false;
	//for(unsigned int c(0);c<N_;c++){
	//for(unsigned int i(0);i<n_;i++){
	//for(unsigned int j(0);j<M_(c)-1;j++){
	//EVec_[c](i,j) = H_(i,j);
	//}
	//}
	//}
	//compute_H();
	//select_eigenvectors(M_(0));
	//for(unsigned int c(0);c<N_;c++){
	//for(unsigned int i(0);i<n_;i++){
	//EVec_[c](i,M_(c)-1) = evec_(i,M_(c));
	//}
	//std::cout<<(EVec_[c].trans_conj()*EVec_[c]).chop()<<std::endl;
	//}
	///*}*/
	/*{Description*/
	compute_H();
	diagonalize_H(H_);
	degenerate_=false;
	for(unsigned int c(0);c<N_;c++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_(c)-1;j++){
				EVec_[c](i,j) = H_(i,j);
			}
		}
	}
	compute_H();
	select_eigenvectors(M_(0));
	Vector<std::complex<double> > v(n_);
	for(unsigned int c(0);c<N_;c++){
		switch(c){
			case 0:
				{
					for(unsigned int i(0);i<n_;i++){
						v(i) = evec_(i,M_(c));
					}
				}break;
			case 1:
				{
					for(unsigned int i(0);i<n_;i++){
						v(i) = evec_(i,M_(c)-1);
					}
				}break;
			case 2:
				{
					for(unsigned int i(0);i<n_;i++){
						v(i) = -evec_(i,M_(c)) - evec_(i,M_(c)-1);
					}
				}break;
		}
		double norm(v.norm());
		for(unsigned int i(0);i<n_;i++){ EVec_[c](i,M_(c)-1) = v(i)/norm; }
		std::cout<<(EVec_[c].trans_conj()*EVec_[c]).chop()<<std::endl;
		std::cout<<"..................."<<std::endl;
		std::cout<<EVec_[c].chop()<<std::endl;
		std::cout<<"..................."<<std::endl;
	}
	/*}*/
}
