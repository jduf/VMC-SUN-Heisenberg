#include "HoneycombSU4.hpp"

HoneycombSU4::HoneycombSU4(Parseur& P):
	Honeycomb<double>(P)
{
	if(!P.status() || N_ != 4){
		if(m_==Ly_*Lx_){
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
				std::string filename("honeycomb");
				filename +="-N" + tostring(N_);
				filename +="-S" + tostring(n_);
				filename += "-" + tostring(Lx_) +"x"+ tostring(Ly_);
				if(bc_ == 1){ filename += "-P";} 
				else { filename += "-A";}
				save(filename);
			} else {
				if(bc_ == 1){
					std::cerr<<"HoneycombSU4 : degeneate for PBC"<<std::endl;
				} else {
					std::cerr<<"HoneycombSU4 : degeneate for APBC"<<std::endl;
				}
			}
		} else {
			std::cerr<<"HoneycombSU4 : the cluster is not a square"<<std::endl;
		}
	}
}

HoneycombSU4::~HoneycombSU4(){}

void HoneycombSU4::compute_T(){
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
	//std::cout<<" | ";
	//for(unsigned int j(0);j<T_.col();j++){
		//std::cout<<j<<" ";
	//}
	//std::cout<<std::endl;
	//for(unsigned int j(0);j<T_.col();j++){
		//std::cout<<"__";
	//}
	//std::cout<<std::endl;
	//for(unsigned int i(0);i<T_.row();i++){
		//std::cout<<i<<"| ";
		//for(unsigned int j(0);j<T_.col();j++){
			//std::cout<<T_(i,j)<<" ";
		//}
		//std::cout<<std::endl;
	//}
	T_ += T_.transpose();
}

void HoneycombSU4::save(std::string filename){
	Write w(filename+".jdbin");
	RST rst;
	rst.title("SU(4) Honeycomb lattice","~");
	rst.item("pi-flux configuration");
	rst.item("4 sites per unit cell");
	rst.item("Lx x Ly unit cells");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("wf (wave function)",wf_);
	w("N (N of SU(N))",N_);
	w("m (number of unit cell)",m_);
	w("sts (connected sites)",sts_);
	w("EVec (unitary matrix)",EVec_);
	w("bc (boundary condition)",bc_);
	w("Lx (x-dimension)",Lx_);
	w("Ly (y-dimension)",Ly_);
}

