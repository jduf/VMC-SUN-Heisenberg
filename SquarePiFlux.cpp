#include "SquarePiFlux.hpp"

SquarePiFlux::SquarePiFlux(Parseur& P):
	Square<std::complex<double> >(P)
{
	if(!P.status()){
		if(n_==Ly_*Lx_){
			compute_EVec();
			for(unsigned int spin(0);spin<N_;spin++){
				for(unsigned int i(0);i<n_;i++){
					for(unsigned int j(0);j<m_;j++){
						EVec_(i+spin*n_,j) = T_(i,j);
					}
				}
			}
			if(successful_){
				std::string filename("square-piflux");
				filename += "-N" + tostring(N_);
				filename += "-S" + tostring(n_);
				filename += "-" + tostring(Ly_) + "x" + tostring(Lx_);
				if(bc_ == 1){ filename += "-P";} 
				else { filename += "-A";}
				save(filename);
			} else {
				std::cerr<<"SquarePiFlux : degeneate"<<std::endl;
			}
		} else {
			std::cerr<<"SquarePiFlux : the cluster is not a SquarePiFlux"<<std::endl;
		}
	}
}

SquarePiFlux::~SquarePiFlux(){}

void SquarePiFlux::compute_EVec(){
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
	std::cerr<<"SquarePiFlux : compute_EVec : modified the flux disposition..."<<std::endl;
	std::cout<<T_<<std::endl;
	T_ += T_.trans_conj(); 
	diagonalize_EVec('H');
}

void SquarePiFlux::save(std::string filename){
	Write w(filename+".jdbin");
	RST rst;
	rst.text("Chiral spin liquid, with 2pi/N flux per plaquette");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("is_complex",true);
	w("N_",N_);
	w("m_",m_);
	w("sts",sts_);
	w("EVec",EVec_);
	w("bc_",bc_);
	w("Ly_",Ly_);
	w("Lx_",Lx_);
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
