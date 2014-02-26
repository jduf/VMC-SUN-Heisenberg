#include "SquarePiFlux.hpp"

SquarePiFlux::SquarePiFlux(unsigned int N, unsigned int n, unsigned int m):
	Square<std::complex<double> >(N,n,m,"square-csl")
{
	rst_.text("Chiral spin liquid, with 2pi/N flux per plaquette");
}

SquarePiFlux::~SquarePiFlux(){}

void SquarePiFlux::compute_T(){
	double t(-1.0);
	double phi(2*M_PI/N_);
	for(unsigned int i(0); i< Ly_; i++){
		for(unsigned int j(0); j< Lx_; j++){
			if(j+1 == Lx_){ T_( i*Lx_ , i*Lx_ + j) = bc_*t; }
			else{ T_(i*Lx_ + j , i*Lx_ + j + 1) = t; }
			if(i+1 == Ly_ ){ T_(j, i*Lx_ + j) = std::polar(bc_*t,-((j%N_)+1)*phi); } 
			//if(i+1 == Ly_ ){ T(j, i*Lx_ + j) = std::polar(bc_*t,-(j%N_)*phi); } 
			else{ T_(i*Lx_ + j, (i+1)*Lx_ + j)= std::polar(t,((j%N_)+1)*phi); }
		}
	}
	std::cerr<<"SquarePiFlux : compute_EVec : new use of polar, check that it is correct"<<std::endl;
	std::cerr<<"                            : modified the flux disposition..."<<std::endl;
	T_ += T_.trans_conj(); 
}

void SquarePiFlux::create(double x){
	compute_T();
	diagonalize_T('H');
	for(unsigned int spin(0);spin<N_;spin++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_;j++){
				EVec_(i+spin*n_,j) = T_(i,j);
			}
		}
	}
}

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
