#include "SquarePiFlux.hpp"

SquarePiFlux::SquarePiFlux(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc):
	System(ref,N,m,n,M,bc),
	Square<std::complex<double> >(1,1,N,"square-csl")
{
	if(status_==1){
		init_fermionic();

		system_info_.text("Chiral spin liquid, with 2pi/N flux per plaquette");
		std::cout<<"check everything"<<std::endl;
	}
}

/*{method needed for running*/
void SquarePiFlux::compute_H(){
	double t(-1.0);
	double phi(2*M_PI/N_);
	for(unsigned int i(0); i< Ly_; i++){
		for(unsigned int j(0); j< Lx_; j++){
			if(j+1 == Lx_){ H_( i*Lx_ , i*Lx_ + j) = bc_*t; }
			else{ H_(i*Lx_ + j , i*Lx_ + j + 1) = t; }
			if(i+1 == Ly_ ){ H_(j, i*Lx_ + j) = std::polar(bc_*t,-((j%N_)+1)*phi); } 
			//if(i+1 == Ly_ ){ T(j, i*Lx_ + j) = std::polar(bc_*t,-(j%N_)*phi); } 
			else{ H_(i*Lx_ + j, (i+1)*Lx_ + j)= std::polar(t,((j%N_)+1)*phi); }
		}
	}
	std::cerr<<"SquarePiFlux : compute_EVec : new use of polar, check that it is correct"<<std::endl;
	std::cerr<<"                            : modified the flux disposition..."<<std::endl;
	H_ += H_.trans_conj(); 
}

void SquarePiFlux::create(){
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);

	compute_H();
	diagonalize_H(H_);
	for(unsigned int c(0);c<N_;c++){
			EVec_[c].set(n_,M_(c));
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = H_(i,j);
				}
		}
	}
}
/*}*/

/*{method needed for checking*/
void SquarePiFlux::check(){
	compute_H();
	std::cout<<H_.chop(1e-6)<<std::endl;
}
/*}*/

//{//csl for Vishvanath (uses majorana representation)
//for(unsigned int i(0); i< Ly_; i++){
//for(unsigned int j(0); j< Lx_; j++){
//if(j+1 == Lx_){// x hopping
//H(i*Lx_ , i*Lx_ + j) = t;
//if(i % 2 == 0){
//T(i*Lx_ , i*Lx_ + j) = bc_*t;
//} else {
//T(i*Lx_ , i*Lx_ + j) = -bc_*t;
//}
//} else {
//H( i*Lx_ + j , i*Lx_ + j + 1) = t; 
//if(i % 2 == 0){
//T( i*Lx_ + j , i*Lx_ + j + 1) = t; 
//} else {
//T( i*Lx_ + j , i*Lx_ + j + 1) = -t; 
//}
//}
//if(i+1 == Ly_ ){// y hopping
//H(j, i*Lx_ + j) = t;
//T(j, i*Lx_ + j) = bc_*t;
//} else{
//H(i*Lx_ + j, (i+1)*Lx_ + j) = t;
//T(i*Lx_ + j, (i+1)*Lx_ + j) = t;
//}
//}
//}
//H += H.transpose();
//T += T.transpose();
//} 
