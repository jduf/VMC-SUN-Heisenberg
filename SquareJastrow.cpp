#include "SquareJastrow.hpp"

SquareJastrow::SquareJastrow(Parseur& P):
	Square<double>(P,"square-Jastrow"),
	nu_(P.get<double>("nu")),
	nn_(n_,z_),
	sl_(n_),
	omega_(N_,N_,1.0)
{
	if(!P.status()){
		//if(P.get<bool>("study")){
		//compute_T();
		//band_structure();
		//} else {
		//compute_T();
		//diagonalize_T('S');
		//EVec_.set(N_*n_,n_,0.0);
		//for(unsigned int i(0);i<N_*n_;i++){
		//for(unsigned int j(0);j<n_;j++){
		//EVec_(i,j) = T_(i,j);
		//}
		//}
		//
		//if(successful_){
		//std::cout<<"prob"<<std::endl;
		//}
		////if(successful_){
		//filename_ += "-N" + tostring(N_);
		//filename_ += "-S" + tostring(n_);
		//filename_ += "-" + tostring(Lx_) + "x" + tostring(Ly_);
		//if(bc_ == 1){ filename_ += "-P";} 
		//else { filename_ += "-A";}
		//filename_ += "-AF+" + tostring(P.get<double>("AF"));
		//
		//save();
		////} else {
		////std::cerr<<"SquareJastrow : degeneate"<<std::endl;
		////}
		//}

		filename_ += "-N" + tostring(N_);
		filename_ += "-S" + tostring(n_);
		filename_ += "-" + tostring(Lx_) + "x" + tostring(Ly_);
		if(bc_ == 1){ filename_ += "-P";} 
		else { filename_ += "-A";}
		filename_ += "-nu+" + tostring(nu_);

		compute_nn();
		compute_sublattice();
		compute_omega();
		save();
	} else {
		std::cerr<<"SquareJastrow : need to provide nu"<<std::endl;
	}
}

SquareJastrow::~SquareJastrow(){}

void SquareJastrow::save(){
	Write w(filename_+".jdbin");
	RST rst;
	rst.text("Staggered magnetic field, Becca's idea to mimic an on site chemical potential");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("wf (wave function)",wf_);
	w("N (N of SU(N))",N_);
	w("m (m=n/N)",m_);
	w("bc (boundary condition)",bc_);
	w("Lx (x-dimension)",Lx_);
	w("Ly (y-dimension)",Ly_);
	w("nu (jastrow coefficient)",nu_);
	w("nn (nearst neighbours)",nn_);
	w("sl (sublattice)",sl_);
	w("omega (omega)",omega_);
	w("sts (connected sites)",sts_);
}

void SquareJastrow::compute_nn(){
	unsigned int k(0);
	for(unsigned int i(0);i<n_;i++){
		k=0;
		for(unsigned int j(0);j<n_;j++){
			if(H_(i,j)){ 
				nn_(i,k) = j;
				k++;
			}
		}
	}
}

void SquareJastrow::compute_sublattice(){
	unsigned int k(0);
	for(unsigned int i(0);i<n_;i++){
		sl_(i) = k % N_;
		if((i+1)%Lx_==0){
			k++;
		}
		k++;
	}
}

void SquareJastrow::compute_omega(){
	if(N_==2){
		omega_(1,1) = -1.0;
	}
	if(N_==3){
		omega_(1,1) = std::polar(1.0,2.0*M_PI/3.0);
		omega_(2,2) = std::polar(1.0,2.0*M_PI/3.0);
		omega_(1,2) = std::polar(1.0,4.0*M_PI/3.0);
		omega_(2,1) = std::polar(1.0,4.0*M_PI/3.0);
	}
}
