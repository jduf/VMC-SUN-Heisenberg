#include "SquareFermi.hpp"

SquareFermi::SquareFermi(Parseur& P):
	Square<double>(P)
{
	if(!P.status()){
		if(n_==Ly_*Lx_){
			compute_T();
			compute_P();
			compute_band_structure();
			//for(unsigned int spin(0);spin<N_;spin++){
				//for(unsigned int i(0);i<n_;i++){
					//for(unsigned int j(0);j<m_;j++){
						//EVec(i+spin*n_,j) = T(i,j);
					//}
				//}
			//}
			if(successful_){
				std::string filename("square-fermi");
				filename += "-N" + tostring(N_);
				filename += "-S" + tostring(n_);
				filename += "-" + tostring(Ly_) + "x" + tostring(Lx_);
				if(bc_ == 1){ filename += "-P";} 
				else { filename += "-A";}
				save(filename);
			} else {
				std::cerr<<"SquareFermi : degeneate"<<std::endl;
			}
		} else {
			std::cerr<<"SquareFermi : the cluster is not a square"<<std::endl;
		}
	}
}

SquareFermi::~SquareFermi(){}

void SquareFermi::compute_T(){
	T_.set(n_,n_,0.0);
	double t(-1.0);
	for(unsigned int i(0); i < n_; i++){
		/*horizontal hopping*/
		if( (i+1) % Lx_ ){ T_(i,i+1) = t;}	
		else { T_(i+1-Lx_,i) = bc_*t;}
		/*vertical hopping*/
		if( i+Lx_ < n_ ){  T_(i,i+Lx_) = t; } 
		else { T_(i-(Ly_-1)*Lx_,i) = bc_*t;}
	}
	T_ += T_.transpose();
}

void SquareFermi::compute_P(){
	for(unsigned int i(0); i < n_; i++){
		/*horizontal hopping*/
		if( (i % Ly_)  < Ly_ - 1 ){ Px_(i,i+1) = 1; }
		else { Px_(i,i+1-Lx_) = bc_; }
		/*vertical hopping*/
		if( i+Lx_ < n_ ){ Py_(i,i+Lx_) = 1; }
		else { Py_(i,i-(Ly_-1)*Lx_) = bc_; }
	}
	//std::cout<<T_*Px_-Px_*T_<<std::endl;
	//std::cout<<T_*Py_-Py_*T_<<std::endl;
}


void SquareFermi::save(std::string filename){
	Write w(filename+".jdbin");
	RST rst;
	rst.text("fermi : all colors experience the same Hamiltonian");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("is_complex",false);
	w("N_",N_);
	w("m_",m_);
	w("sts",sts_);
	w("EVec",EVec_);
	w("bc_",bc_);
	w("Ly_",Ly_);
	w("Lx_",Lx_);
}

