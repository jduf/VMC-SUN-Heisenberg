#include "HoneycombSU4.hpp"

HoneycombSU4::HoneycombSU4(unsigned int N, unsigned int n, unsigned int m):
	Honeycomb<double>(N,n,m,"honeycomb-pi-flux")
{
	rst_.title("SU(4) Honeycomb lattice","~");
	rst_.item("pi-flux configuration");
	rst_.item("4 sites per unit cell");
	rst_.item("Lx x Ly unit cells");

	if(bc_ == 1){ filename_ += "-P";} 
	else { filename_ += "-A";}
	if(N_ != 4 ||M_==Ly_*Lx_){
		std::cerr<<"HoneycombSU4 : the cluster is not a square"<<std::endl;
		std::cerr<<"HoneycombSU4 : or N!=4"<<std::endl;
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

void HoneycombSU4::create(double x){
	compute_T();
	diagonalize_T('S');
	for(unsigned int spin(0);spin<N_;spin++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_;j++){
				EVec_(i+spin*n_,j) = T_(i,j);
			}
		}
	}
}

void HoneycombSU4::check(){
	compute_T();
	std::cout<<T_.chop(1e-6)<<std::endl;
}
