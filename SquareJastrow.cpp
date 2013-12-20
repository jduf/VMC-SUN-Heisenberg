#include "SquareJastrow.hpp"

SquareJastrow::SquareJastrow(Parseur& P):
	Square<double>(P,"square-Jastrow"),
	nu_(P.get<double>("nu")),
	nn_(n_,3*z_),
	sl_(n_),
	omega_(N_,N_,1.0)
{
	if(!P.status()){
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

Vector<unsigned int> SquareJastrow::get_neighbourg(unsigned int i){
	Vector<unsigned int> neighbourg(z_);
	/*+x neighbour*/
	if((i+1)%Lx_!=0){ neighbourg(0) = i+1; } 
	else { neighbourg(0) = (i/Lx_)*Lx_; }
	/*+y neighbour*/
	if(i<n_-Lx_){ neighbourg(1) = i+Lx_; }
	else { neighbourg(1) = i-n_+Lx_; }
	/*-x neighbour*/
	if(i%Lx_){ neighbourg(2) = i-1; }
	else { neighbourg(2) = i+Lx_-1; }
	/*-y neighbour*/
	if(i>=Lx_){ neighbourg(3) = i-Lx_; }
	else { neighbourg(3) = n_-Lx_+i; }

	return neighbourg;
}

void SquareJastrow::compute_nn(){
	Vector<unsigned int> neighbourg;
	for(unsigned int i(0);i<n_;i++){
		neighbourg = get_neighbourg(i);
		for(unsigned int j(0);j<z_;j++){
			nn_(i,j) = neighbourg(j);
		}
		unsigned int l(z_);
		for(unsigned int j(0);j<z_;j++){
			neighbourg = get_neighbourg(nn_(i,j));
			for(unsigned int k(j);k<j+2;k++){
				nn_(i,l) = neighbourg(k%z_);
				l++;
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
