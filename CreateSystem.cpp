#include "CreateSystem.hpp"

template<>
void CreateState<double>::compute_H(){
	is_complex = false;
	mat_type = 'S';
	switch(N_n){
		case 2:
			{
				if(N_m % 2 == 0){ H(0,N_site -1 ) = 1.0; }
				else { H(0,N_site -1 ) = -1.0;}
				for(unsigned int i(0); i< N_site-1; i++){
					H(i,i+1) = -1.0;
				}
				H += H.transpose();
				T = H;
				break;
			}
		case 4:
			{
				unsigned int N_row(0),N_col(0);
				N_row = sqrt(N_site);
				N_col = N_site/N_row;
				assert(N_site==N_row*N_col);
				unsigned int j(0);
				for(unsigned int i(0); i< N_site; i++){
					if((i+1) % N_col == 0){H(j*N_col,i) = 1.0;j++;}
					else {H(i,i+1) = -1.0;}
					if(i+N_col < N_site){ H(i,i+N_col) = -1.0;}
					else{ H(i+N_col-N_site,i) = 1.0;}
				}
				H += H.transpose();
				T = H;
				break;
			}
		default:
			{
				std::cerr<<"init_H<real> : lattice type undefined"<<std::endl;
			}
	}
}

template<>
void CreateState<std::complex<double> >::compute_H(){
	is_complex = true;
	mat_type = 'H';
	unsigned int N_site(H.size());
	switch(N_n){
		case 4:
			{
				unsigned int N_row(0),N_col(0);
				N_row = sqrt(N_site);
				N_col = N_site/N_row;
				assert(N_site==N_row*N_col);
				double phi(2*M_PI/N_spin);
				unsigned int j(0);
				for(unsigned int i(0); i< N_site; i++){
					if((i+1) % N_col == 0){
						T(j*N_col,i) = 1.0;
						H(j*N_col,i) = 1.0;
						j++;
					} else {
						T(i,i+1) = -1.0;
						H(i,i+1) = -1.0;
					}
					if(i+N_col < N_site){
						T(i,i+N_col) = -std::polar(1.0,((i%N_spin)+1)*phi);
						H(i,i+N_col) = -1.0;
					} else{
						T(i+N_col-N_site,i) = std::polar(1.0,-((i%N_spin)+1)*phi);
						H(i+N_col-N_site,i) = 1.0; 
					}
				}
				H += H.transpose();
				T += T.trans_conj(); 
				break;
			}
		default:
			{
				std::cerr<<"init_H<complex> : lattice type undefined"<<std::endl;
			}
	}
}

