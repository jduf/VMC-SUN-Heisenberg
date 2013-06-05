#include "Square.hpp"

Square::Square(Parseur& P):
	CreateSystem<std::complex<double> >(P,4),
	N_row(floor(sqrt(N_site))),
	N_col(floor(sqrt(N_site)))
{
	P.set("bc",bc);
	if(N_site==N_row*N_col){
		is_complex=true;
		mat_type='H';
		compute_T();
		compute_sts();
		compute_EVec();
		if(successful){
			filename = "square-csl-N"+tostring(N_spin) + "-S" + tostring(N_site);
			if(bc == 1){
				std::cout<<"system with PBC created"<<std::endl;
				filename += "-" + tostring(N_row) +"x"+ tostring(N_col) + "-P" ;
			} else {
				std::cout<<"system with APBC created"<<std::endl;
				filename += "-" + tostring(N_row) +"x"+ tostring(N_col) + "-A" ;
			}
			rst.text("Chiral spin liquid, with 2pi/N flux per plaquette");
			rst.np();
			save();
			(*w)("N_row",N_row);
			(*w)("N_col",N_col);
		} else {
			if(bc == 1){
				std::cerr<<"CreateSystem : degeneate for PBC"<<std::endl;
			} else {
				std::cerr<<"CreateSystem : degeneate for APBC"<<std::endl;
			}
		}
	} else {
		std::cerr<<"CreateSystem : the cluster is not a square"<<std::endl;
	}
}

Square::~Square()
{}

void Square::compute_T(){
	double t(-1.0);
	double phi(2*M_PI/N_spin);
	for(unsigned int i(0); i< N_row; i++){
		for(unsigned int j(0); j< N_col; j++){
			if(j+1 == N_col){
				H( i*N_col , i*N_col + j) = t;
				T( i*N_col , i*N_col + j) = bc*t;
			} else { 
				H( i*N_col + j , i*N_col + j + 1) = t; 
				T( i*N_col + j , i*N_col + j + 1) = t; 
			}
			if(i+1 == N_row ){
				H(j, i*N_col + j) = t;
				T(j, i*N_col + j) = bc*t*std::polar(1.0,-((j%N_spin)+1)*phi);
			} else{
				H(i*N_col + j, (i+1)*N_col + j) = t;
				T(i*N_col + j, (i+1)*N_col + j) = t*std::polar(1.0,((j%N_spin)+1)*phi);
			}
		}
	}
	H += H.transpose();
	T += T.trans_conj(); 
}

	//{//csl for Vishvanath (uses majorana representation)
		//for(unsigned int i(0); i< N_row; i++){
			//for(unsigned int j(0); j< N_col; j++){
				//if(j+1 == N_col){// x hopping
					//H(i*N_col , i*N_col + j) = t;
					//if(i % 2 == 0){
						//T(i*N_col , i*N_col + j) = bc*t;
					//} else {
						//T(i*N_col , i*N_col + j) = -bc*t;
					//}
				//} else {
					//H( i*N_col + j , i*N_col + j + 1) = t; 
					//if(i % 2 == 0){
						//T( i*N_col + j , i*N_col + j + 1) = t; 
					//} else {
						//T( i*N_col + j , i*N_col + j + 1) = -t; 
					//}
				//}
				//if(i+1 == N_row ){// y hopping
					//H(j, i*N_col + j) = t;
					//T(j, i*N_col + j) = bc*t;
				//} else{
					//H(i*N_col + j, (i+1)*N_col + j) = t;
					//T(i*N_col + j, (i+1)*N_col + j) = t;
				//}
			//}
		//}
		//H += H.transpose();
		//T += T.transpose();
	//} 
