#include "SquarePiFlux.hpp"

SquarePiFlux::SquarePiFlux(Parseur& P):
	Square<std::complex<double> >(P,"square-csl")
{
	if(!P.status()){
		if(n_==Ly_*Lx_){
			compute_T();

			diagonalize_T('H');
			for(unsigned int spin(0);spin<N_;spin++){
				for(unsigned int i(0);i<n_;i++){
					for(unsigned int j(0);j<m_;j++){
						EVec_(i+spin*n_,j) = T_(i,j);
					}
				}
			}
			if(successful_){
				filename_ += "-N" + tostring(N_);
				filename_ += "-S" + tostring(n_);
				filename_ += "-" + tostring(Lx_) + "x" + tostring(Ly_);
				if(bc_ == 1){ filename_ += "-P";} 
				else { filename_ += "-A";}
				save();
			} else {
				std::cerr<<"SquarePiFlux : degeneate"<<std::endl;
			}
		} else {
			std::cerr<<"SquarePiFlux : the cluster is not a SquarePiFlux"<<std::endl;
		}
	}
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
	//std::cout<<T_<<std::endl;
	T_ += T_.trans_conj(); 
}

void SquarePiFlux::save(){
	Write w(filename_+".jdbin");
	RST rst;
	rst.text("Chiral spin liquid, with 2pi/N flux per plaquette");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("wf (wave function)",wf_);
	w("N (N of SU(N))",N_);
	w("m (m=n/N)",m_);
	w("bc (boundary condition)",bc_);
	w("Lx (x-dimension)",Lx_);
	w("Ly (y-dimension)",Ly_);
	w("sts (connected sites)",sts_);
	w("EVec (unitary matrix)",EVec_);
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
