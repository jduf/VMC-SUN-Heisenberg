#include "CreateSystem.hpp"

template<>
void CreateSystem<double>::compute_H(unsigned int N_row, unsigned int N_col, double parity){
	std::cout<<N_site<<" "<<N_col<<" "<<N_row<<" "<<parity<<std::endl;
	H.set(0.0);
	T.set(0.0);
	is_complex = false;
	mat_type = 'S';
	for(unsigned int i(0); i< N_row; i++){
		for(unsigned int j(0); j< N_col; j++){
			if(j+1 == N_col){ H(i*N_col , i*N_col + j) = parity; }
			else { H( i*N_col + j , i*N_col + j + 1) = -1.0; }
			if(i+1 == N_row ){ H(j, i*N_col + j) = parity;}
			else{ H(i*N_col + j, (i+1)*N_col + j) = -1.0;}
		}
	}
	H += H.transpose();
	T = H;
}

template<>
void CreateSystem<std::complex<double> >::compute_H(unsigned int N_row, unsigned int N_col, double parity){
	std::cout<<N_site<<" "<<N_col<<" "<<N_row<<" "<<parity<<std::endl;
	H.set(0.0);
	T.set(0.0);
	is_complex = true;
	mat_type = 'H';
	unsigned int N_site(H.size());
	double phi(2*M_PI/N_spin);
	unsigned int j(0);
	for(unsigned int i(0); i< N_row; i++){
		for(unsigned int j(0); j< N_col; j++){
			if(j+1 == N_col){
				H( i*N_col , i*N_col + j) = parity;
				T( i*N_col , i*N_col + j) = parity;
			} else { 
				H( i*N_col + j , i*N_col + j + 1) = -1.0; 
				T( i*N_col + j , i*N_col + j + 1) = -1.0; 
			}
			if(i+1 == N_row ){
				H(j, i*N_col + j) = parity;
				T(j, i*N_col + j) = parity*std::polar(1.0,-((j%N_spin)+1)*phi);
			} else{
				H(i*N_col + j, (i+1)*N_col + j) = -1.0;
				T(i*N_col + j, (i+1)*N_col + j) = -std::polar(1.0,((j%N_spin)+1)*phi);
			}
		}
	}
	H += H.transpose();
	T += T.trans_conj(); 
	std::cout<<T<<std::endl;
}

template<>
void CreateSystem<double>::compute_H(){
	is_complex = false;
	mat_type = 'S';
	if(N_m % 2 == 0){ 
		filename += "-A";
		H(0,N_site -1 ) = 1.0;
	} else {
		filename += "-P";
		H(0, N_site -1 ) = -1.0;
	}
	for(unsigned int i(0); i< N_site-1; i++){
		H(i,i+1) = -1.0;
	}
	H += H.transpose();
	T = H;
}

template<>
void CreateSystem<std::complex<double> >::compute_H(){
	std::cerr<<"CreateSystem : compute_H not defined for a complex without specifing N_row and N_col"<<std::endl;
}
