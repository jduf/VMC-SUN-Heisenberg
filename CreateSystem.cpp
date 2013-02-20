#include "CreateSystem.hpp"

template<>
void CreateState<double>::init(){
	mat_type = 'S';
	switch(N_n){
		case 2:
			{
				H(0,N_site -1 ) = 1.0;
				for(unsigned int i(0); i< N_site-1; i++){
					H(i,i+1) = -1.0;
				}
				H += H.transpose();
				break;
			}
		case 4:
			{
				unsigned int N_row(0),N_col(0);
				std::cout<<"N_site="<<N_site<<std::endl;
				while(N_row*N_col != N_site || N_row<3 || N_col<3){
					std::cout<<"N_col=";
					std::cin>>N_col;
					std::cout<<"N_row=";
					std::cin>>N_row;
				}
				unsigned int j(0);
				for(unsigned int i(0); i< N_site; i++){
					if((i+1) % N_col == 0){H(j*N_col,i) = 1.0;j++;}
					else {H(i,i+1) = -1.0;}
					if(i+N_col < N_site){ H(i,i+N_col) = -1.0;}
					else{ H(i+N_col-N_site,i) = 1.0;}
				}
				H += H.transpose();
				break;
			}
		default:
			{
				std::cerr<<"init_H<real> : lattice type undefined"<<std::endl;
			}
	}
	filename += 
	w.open(filename);
	w<<false<<N_spin<<N_m<<N_n;
}

template<>
void CreateState<std::complex<double> >::init(){
	mat_type = 'H';
	w<<true<<N_spin<<N_m<<N_n;
	unsigned int N_site(H.size());
	switch(N_n){
		case 4:
			{
				unsigned int N_row(4),N_col(6);
				assert(N_row*N_col==N_site);
				double phi(2*M_PI/N_col);
				unsigned int j(0);
				for(unsigned int i(0); i< N_site; i++){
					if((i+1) % N_col == 0){ H(j*N_col,i) = 1.0; j++;}
					else { H(i,i+1) = -1.0;}
					if(i+N_col < N_site){ H(i,i+N_col) = std::polar(1.0,-((i%N_col)+1)*phi);}
					else{ H(i+N_col-N_site,i) = std::polar(1.0,((i%N_col)+1)*phi);}
				}
				H += H.trans_conj();
				break;
			}
		default:
			{
				std::cerr<<"init_H<complex> : lattice type undefined"<<std::endl;
			}
	}
}

