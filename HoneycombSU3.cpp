#include "HoneycombSU3.hpp"

HoneycombSU3::HoneycombSU3(Parseur& P):
	Honeycomb<double>(P,"honeycomb")
{
	if(!P.status() || N_ != 3){
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
				filename_ +="-N" + tostring(N_);
				filename_ +="-S" + tostring(n_);
				filename_ += "-" + tostring(Lx_) +"x"+ tostring(Ly_);
				if(bc_ == 1){ filename_ += "-P";} 
				else { filename_ += "-A";}
				save();
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
	w("wf (wave function)",wf_);
	w("N (N of SU(N))",N_);
	w("m (number of unit cell)",m_);
	w("bc (boundary condition)",bc_);
	w("Lx (x-dimension)",Lx_);
	w("Ly (y-dimension)",Ly_);
	w("sts (connected sites)",sts_);
	w("EVec (unitary matrix)",EVec_);
}

//{//for SU(3) with 12 sites p00-flux
//unsigned int i(0);
//unsigned int j(0);
//for(unsigned int ix(0);ix<Lx_;ix++){
//i = j;
//j = Ly_*(2*ix+1)*6;
//for(unsigned int iy(0);iy<Ly_;iy++){ 
////cell (A)
//H_(i,j) = t;
//H_(i,i+1) = t;
//H_(j,j+1) = t;
//T_(i,j) = th;
//T_(i,i+1) = th;
//T_(j,j+1) = th;
//i++;
//j++;
//if( j+Ly_*6>n_ ){
//H_(j%(Ly_*6),j)= bc_*t; 
//T_(j%(Ly_*6),j)= bc_*td; 
//} else { 
//H_(j,j+Ly_*6) = t; 
//T_(j,j+Ly_*6) = td; 
//}
//H_(i,i+1) = t;
//H_(j,j+1) = t;
//T_(i,i+1) = th;
//T_(j,j+1) = th;
//i++;
//j++;			
////cell (B)
//H_(i,j) = t;
//H_(i,i+1) = t;
//H_(j,j+1) = t;
//T_(i,j) = th;
//T_(i,i+1) = td;
//T_(j,j+1) = td;
//i++;
//j++;
//if( j+Ly_*6>n_ ){
//H_(j%(Ly_*6),j)= bc_*t; 
//T_(j%(Ly_*6),j)= bc_*th; 
//} else {
//H_(j,j+Ly_*6) = t; 
//T_(j,j+Ly_*6) = th; 
//}
//H_(i,i+1) = t;
//H_(j,j+1) = t;
//T_(i,i+1) = th;
//T_(j,j+1) = th;
//i++;
//j++;			
////cell (C)
//H_(i,j) = t;
//H_(i,i+1) = t;
//H_(j,j+1) = t;
//T_(i,j) = td;
//T_(i,i+1) = th;
//T_(j,j+1) = th;
//i++;
//j++;
//if( j+Ly_*6>n_ ){
//H_(j%(Ly_*6),j)= bc_*t; 
//T_(j%(Ly_*6),j)= bc_*th; 
//} else{
//H_(j,j+Ly_*6) = t;
//T_(j,j+Ly_*6) = th;
//}
//if( iy+1 == Ly_  ){ // link (A)-(C)
//H_(i+1,j)= bc_*t; 
//H_((i+1)-Ly_*6,i) = bc_*t; 
//T_(i+1,j)= bc_*td; 
//T_((i+1)-Ly_*6,i) = bc_*td; 
//} else { // link (A)-(C)
//H_(i,i+1) = t; 
//H_(j,j+1) = t; 
//T_(i,i+1) = td; 
//T_(j,j+1) = td; 
//}
//i++;
//j++;
//}
//}
//}
//{//for SU(3) (0pipi) with 6 sites
//unsigned int i(0);
//for(unsigned int l(0);l<Ly_;l++){
//for(unsigned int c(0);c<Lx_;c++){
////std::cout<<i<<std::endl;
//H_(i,i+1) = t; //0-1
//T_(i,i+1) = th; //0-1
//i++;//1
//H_(i,i+1) = t; //1-2
//T_(i,i+1) = th; //1-2
//i++;//2
//H_(i,i+1) = t; //2-3
//T_(i,i+1) = td; //2-3
//H_(i,i+2) = t; //2-4
//T_(i,i+2) = th; //2-4
//i++;//3
//H_(i,i+2) = t; //3-5
//T_(i,i+2) = th; //3-5
//if(c+1==Lx_){
//H_((i+3) - Lx_*6, i) = t; // 0-3 b
//T_((i+3) - Lx_*6, i) = bc_*th; // 0-3 b
//} else {
//H_(i,i+3) = t; //3-6
//T_(i,i+3) = th; //3-6
//}
//i++;//4
//if(l+1==Ly_){
//H_(i-4-l*Lx_*6,i) = t; // 0-4 b
//T_(i-4-l*Lx_*6,i) = bc_*td; // 0-4 b
//} else {
//H_(i,i-4+Lx_*6) = t; // 4-18
//T_(i,i-4+Lx_*6) = td; // 4-18
//}
//i++;//5
//if(l+1==Ly_){
//H_(i-4-l*Lx_*6,i) = t; // 1-5 b
//T_(i-4-l*Lx_*6,i) = bc_*td; // 1-5 b
//} else {
//H_(i,i-4+Lx_*6) = t; // 5-19
//T_(i,i-4+Lx_*6) = td; // 5-19
//}
////std::cout<<i<<std::endl;
//if(c+1==Lx_){
//H_((i+5) - Lx_*6, i) = t; // 4-5 b
//T_((i+5) - Lx_*6, i) = bc_*th; // 4-5 b
//} else {
//H_(i,i+5) = t; //5-10
//T_(i,i+5) = th; //5-10
//}
//i++;//6
//}
//}
//}
//{//for SU(3) si cellule pi au milieu (vérifier que tout fonctionne) with 6 sites
//unsigned int i(0);
//for(unsigned int l(0);l<Ly_;l++){
//for(unsigned int c(0);c<Lx_;c++){
////std::cout<<i<<std::endl;
////0
//H_(i,i+1) = t; //0-1
//T_(i,i+1) = td; 
//i++;//1
//H_(i,i+1) = t; //1-2
//T_(i,i+1) = th; 
//if(l != 0 && l+1 != Ly_ && c+1 != Lx_){
//H_(i-9,i) = t; //1-4
//T_(i-9,i) = th; 
//} else {
//if(c+1==Lx_){
//if(l==0){
//H_(i,n_-2-c*6) = t; // 1-4 bc_ bc_
//T_(i,n_-2-c*6) = th; // 1-4 bc_ bc_
//} else {
//H_(i+3-(2*Lx_-1)*6,i) = t; //1-4
//T_(i+3-(2*Lx_-1)*6,i) = th; 
//}
//} else {
//if(l==0){
//H_(i,i+3+((Ly_-1)*Lx_+1)*6) = t; //1-4 bc_
//T_(i,i+3+((Ly_-1)*Lx_+1)*6) = bc_*th; 
//} else {//l+1==Ly_
//H_(i-9,i) = t; //1-4
//T_(i-9,i) = th; 
//}
//
//}
//}
//i++;//2
////std::cout<<i<<std::endl;
//H_(i,i+1) = t; //2-3
//T_(i,i+1) = td;
//if(c+1 == Lx_){
//H_(i+3-c*6,i) = t; //2-5 bc_
//T_(i+3-c*6,i) = bc_*th;
//} else {
//H_(i,i+9) = t;//2-5
//T_(i,i+9) = th;
//}
//i++;//3
//H_(i,i+1) = t; //3-4
//T_(i,i+1) = th;
//if(l+1==Ly_){
//H_(i-3-l*Lx_*6,i) = t;//3-0 bc_
//T_(i-3-l*Lx_*6,i) = bc_*th;
//} else {
//H_(i,i-3+Lx_*6) = t;//3-0
//T_(i,i-3+Lx_*6) = th;
//}
//i++;//4
//H_(i,i+1) = t; //4-5
//T_(i,i+1) = td; 
//i++;//5
//H_(i-5,i) = t; //0-5
//T_(i-5,i) = th; 
////std::cout<<i<<std::endl;
//i++;//6
//}
//}
//}
//{
//for(unsigned int l(0);l<Ly_;l++){
//for(unsigned int c(0);c<Lx_;c++){
//std::cout<<i<<std::endl;
//H_(i,i+1) = t; //0-1
//T_(i,i+1) = th; //0-1
//i++;//1
//H_(i,i+1) = t; //1-2
//T_(i,i+1) = th; //1-2
//i++;//2
//H_(i,i+1) = t; //2-3
//T_(i,i+1) = td; //2-3
//H_(i,i+2) = t; //2-4
//T_(i,i+2) = th; //2-4
//i++;//3
//H_(i,i+2) = t; //3-5
//T_(i,i+2) = th; //3-5
//if(c+1==Lx_){
//H_((i+3) - Lx_*6, i) = t; // 0-3 b
//T_((i+3) - Lx_*6, i) = bc_*th; // 0-3 b
//} else {
//H_(i,i+3) = t; //3-6
//T_(i,i+3) = th; //3-6
//}
//i++;//4
//if(l+1==Ly_){
//H_(i-4-l*Lx_*6,i) = t; // 0-4 b
//T_(i-4-l*Lx_*6,i) = bc_*td; // 0-4 b
//} else {
//H_(i,i-4+Lx_*6) = t; // 4-18
//T_(i,i-4+Lx_*6) = td; // 4-18
//}
//i++;//5
//if(l+1==Ly_){
//H_(i-4-l*Lx_*6,i) = t; // 1-5 b
//T_(i-4-l*Lx_*6,i) = bc_*td; // 1-5 b
//} else {
//H_(i,i-4+Lx_*6) = t; // 5-19
//T_(i,i-4+Lx_*6) = td; // 5-19
//}
//std::cout<<i<<std::endl;
//if(c+1==Lx_){
//H_((i+5) - Lx_*6, i) = t; // 4-5 b
//T_((i+5) - Lx_*6, i) = bc_*th; // 4-5 b
//} else {
//H_(i,i+5) = t; //5-10
//T_(i,i+5) = th; //5-10
//}
//i++;//6
//}
//}
//}
