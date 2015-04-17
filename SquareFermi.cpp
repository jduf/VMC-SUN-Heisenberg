#include "SquareFermi.hpp"

template<>
void SquareFermi<double>::create(){
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);

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
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);

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

//if(degenerate_){
//Matrix<double> tmp(H_);
//compute_H();
//H_(1,0) -= 0.0001;
//H_(0,1) -= 0.0001;
//Vector<double> eval0;
//Lapack<double>(H_,false,'S').eigensystem(eval0,true);
//unsigned int c(0);
//unsigned int a(this->M_(c)-1);
//unsigned int b(a);
//do{b++;} while (b+1<this->n_ && my::are_equal(eval0(b),eval0(b-1)));
//if(b!=this->M_(c)){ while(a>0 && my::are_equal(eval0(a-1),eval0(a))){a--;} }
//for(unsigned int c(0);c<N_;c++){
//for(unsigned int i(0);i<n_;i++){
//for(unsigned int j(a);j<M_(c);j++){
//EVec_[c](i,j) = H_(i,j);
//}
//}
//}
//std::cout<<H_.chop()<<std::endl;
//std::cout<<std::endl;
//std::cout<<tmp.chop()<<std::endl;
//std::cout<<std::endl;
//std::cout<<(H_-tmp).chop()<<std::endl;
//degenerate_ = false;
//}
