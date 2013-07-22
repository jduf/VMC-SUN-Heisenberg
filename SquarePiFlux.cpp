#include "SquarePiFlux.hpp"

SquarePiFlux::SquarePiFlux(Parseur& P):
	Square<std::complex<double> >(P)
{
	if(!P.status()){
		if(N_site==N_row*N_col){
			mat_type='H';
			compute_T();
			compute_EVec();
			for(unsigned int spin(0);spin<N_spin;spin++){
				for(unsigned int i(0);i<N_site;i++){
					for(unsigned int j(0);j<N_m;j++){
						EVec(i+spin*N_site,j) = T(i,j);
					}
				}
			}
			if(successful){
				std::string filename("square-piflux");
				filename += "-N" + tostring(N_spin);
				filename += "-S" + tostring(N_site);
				filename += "-" + tostring(N_row) + "x" + tostring(N_col);
				if(bc == 1){ filename += "-P";} 
				else { filename += "-A";}
				save(filename);
			} else {
				std::cerr<<"CreateSystem : degeneate"<<std::endl;
			}
		} else {
			std::cerr<<"CreateSystem : the cluster is not a SquarePiFlux"<<std::endl;
		}
	}
}

SquarePiFlux::~SquarePiFlux(){}

void SquarePiFlux::compute_T(){
	double t(-1.0);
	double phi(2*M_PI/N_spin);
	for(unsigned int i(0); i< N_row; i++){
		for(unsigned int j(0); j< N_col; j++){
			if(j+1 == N_col){ T( i*N_col , i*N_col + j) = bc*t; }
			else { T( i*N_col + j , i*N_col + j + 1) = t; }
			if(i+1 == N_row ){ T(j, i*N_col + j) = bc*t*std::polar(1.0,-((j%N_spin)+1)*phi); } 
			else{ T(i*N_col + j, (i+1)*N_col + j) = t*std::polar(1.0,((j%N_spin)+1)*phi); }
		}
	}
	T += T.trans_conj(); 
}

void SquarePiFlux::save(std::string filename){
	Write w(filename+".jdbin");
	RST rst;
	rst.text("Chiral spin liquid, with 2pi/N flux per plaquette");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("is_complex",true);
	w("N_spin",N_spin);
	w("N_m",N_m);
	w("sts",sts);
	w("EVec",EVec);
	w("bc",bc);
	w("N_row",N_row);
	w("N_col",N_col);
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
