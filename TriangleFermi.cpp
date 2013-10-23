#include "TriangleFermi.hpp"

TriangleFermi::TriangleFermi(Parseur& P):
	Triangle<double>(P)
{
	if(!P.status()){
		compute_T();
		diagonalize_T('S');
		for(unsigned int spin(0);spin<N_;spin++){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<m_;j++){
					EVec_(i+spin*n_,j) = T_(i,j);
				}
			}
		}
		if(successful_){
			std::string filename("triangle-fermi");
			filename += "-N" + tostring(N_);
			filename += "-S" + tostring(n_);
			filename += "-" + tostring(Lx_) + "x" + tostring(Ly_);
			if(bc_ == 1){ filename += "-P";} 
			else { filename += "-A";}
			save(filename);
		} else {
			std::cerr<<"TriangleFermi : degeneate"<<std::endl;
		}
	}
}

TriangleFermi::~TriangleFermi(){}

void TriangleFermi::compute_T(){
	double t(-1.0);
	for(unsigned int i(0); i < n_; i++){
		/*horizontal hopping*/
		if( (i+1) % Lx_ ){ T_(i,i+1) = t;}	
		else { T_(i+1-Lx_,i) = bc_*t;}
		/*vertical hopping*/
		if( i+Lx_ < n_ ){  T_(i,i+Lx_) = t; } 
		else { T_(i-(Ly_-1)*Lx_,i) = bc_*t;}
		/*diagonal hopping*/
		if( (i+1) % Lx_ && i+Lx_ < n_ ){  T_(i,i+Lx_+1) = t; } 
		else {
			if(i+1 < n_ ){
				if( !((i+1) % Lx_) ){ T_(i,i+1) = bc_*t;}/*x jump across boundary*/
				if( i+Lx_ >= n_ ){  T_(i-Lx_*(Ly_-1)+1,i) = bc_*t; }/*y jump across boundary*/
			} else {
				T_(0,n_-1) = bc_*bc_*t;
			}
		}
	}
	T_ += T_.transpose();
}

void TriangleFermi::save(std::string filename){
	Write w(filename+".jdbin");
	RST rst;
	rst.text("fermi : all colors experience the same Hamiltonian");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("wf (wave function)",wf_);
	w("N (N of SU(N))",N_);
	w("m (number of unit cell)",m_);
	w("bc (boundary condition)",bc_);
	w("Lx (x-dimension)",Lx_);
	w("Ly (y-dimension)",Ly_);
	w("sts (connected sites)",sts_);
	w("EVec (unitary matrix)",EVec_);
}
