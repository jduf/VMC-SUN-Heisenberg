#include "HoneycombSU3.hpp"

HoneycombSU3::HoneycombSU3(Parseur& P):
	Honeycomb<double>(P,"honeycomb")
{
	ref_(1) = 1;
	ref_(2) = 1;
	if(bc_ == 1){ filename_ += "-P";} 
	else { filename_ += "-A";}
	if(!P.status() || N_ != 3){
		if(M_==Ly_*Lx_){
			compute_T();

			diagonalize_T('S');
			for(unsigned int spin(0);spin<N_;spin++){
				for(unsigned int i(0);i<n_;i++){
					for(unsigned int j(0);j<M_;j++){
						EVec_(i+spin*n_,j) = T_(i,j);
					}
				}
			}
		} else {
			std::cerr<<"HoneycombSU4 : the cluster is not a square"<<std::endl;
		}
	}
}

HoneycombSU3::~HoneycombSU3(){}

void HoneycombSU3::compute_T(){
	double th(1.0);
	double td(-1.0);
	unsigned int i(0);
	for(unsigned int l(0);l<Ly_;l++){
		for(unsigned int c(0);c<Lx_;c++){
			//0
			T_(i,i+1) = td;
			if(l+1<Ly_){
				T_(i,i+1+Lx_*4) = th;
			} else {
				T_(i+1-l*Lx_*4,i) = th*bc_;
			}
			if(c==0){
				T_(i,i+Lx_*4-1) = th*bc_;
			} else {
				T_(i-1,i) = th;
			}
			i+=2;//2
			T_(i,i-1) = th;
			T_(i,i+1) = th; 
			if(l==0){
				T_(i,i+1+(Ly_-1)*Lx_*4) = th*bc_;
			} else {
				T_(i,i+1-Lx_*4) = th;
			}
			i+=2;//4
		}
	}
	T_ += T_.transpose();
}

void HoneycombSU3::save(){
	Write w(filename_+".jdbin");
	std::cerr<<"detail what kind of honeycomb it is"<<std::endl;
	RST rst;
	rst.text("HoneycombSU4 ");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("ref (wave function)",ref_);
	w("n (particles' number)",n_);
	w("N (N of SU(N))",N_);
	w("m (particles per site' number)",m_);
	w("M (particles' number of each color)",M_);
	w("sts (connected sites)",sts_);
	w("EVec (unitary matrix)",EVec_);
	w("bc (boundary condition)",bc_);
	w("Lx (x-dimension)",Lx_);
	w("Ly (y-dimension)",Ly_);
}

