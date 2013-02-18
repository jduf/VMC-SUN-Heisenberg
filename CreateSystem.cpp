#include "CreateSystem.hpp"


void init_H(Matrice<double>& H, unsigned int N_n){
	unsigned int N_site(H.size());
	switch(N_n){
		case 2:
			{
				H(0,1)=-1.0;
				H(0,N_site-1)=1.0;
				for(unsigned int i(1); i< N_site-1; i++){
					H(i,i-1) = -1.0;
					H(i,i+1) = -1.0;
				}
				H(N_site-1,0)=1.0;
				H(N_site-1,N_site-2)=-1.0;
				break;
			}
		case 4:
			{
				unsigned int Nx(3),Ny(4);
				unsigned int j(0);
				for(unsigned int i(0); i< N_site-1; i++){
					if((i+1) % Nx == 0){H(j*Nx,i) = 1.0;j++;}
					else {H(i,i+1) = -1.0;}
					if(i+Nx < N_site){ H(i,i+Nx) = -1.0;}
					else{ H(i+Nx-N_site,i) = 1.0;}
				}
				H(Nx-1,N_site-1) = 1.0;
				H(N_site-Nx,N_site-1) = 1.0;
				H = H+H.transpose();
				break;
			}
		default:
			{
				std::cerr<<"init_H<real> : lattice type undefined"<<std::endl;
			}
	}
}

void init_H(Matrice<std::complex<double> >& H, unsigned int N_n){
	unsigned int N_site(H.size());
	switch(N_n){
		case 4:
			{
				unsigned int Nx(3),Ny(4);
				double phi(M_PI/3.0);
				unsigned int j(0);
				for(unsigned int i(0); i< N_site-1; i++){
					if((i+1) % Nx == 0){ H(j*Nx,i) = 1.0; j++;}
					else { H(i,i+1) = -1.0;}
					if(i+Nx < N_site){ H(i,i+Nx) = -std::polar(1.0,((i%Nx)+1)*phi);}
					else{ H(i+Nx-N_site,i) =std::polar(1.0,((i%Nx)+1)*phi);}
				}
				H(Nx-1,N_site-1) = std::polar(1.0,(((N_site-1)%Nx)+1)*phi);
				H(N_site-Nx,N_site-1) = 1.0;
				H = H+H.trans_conj();
				break;
			}
		default:
			{
				std::cerr<<"init_H<complex> : lattice type undefined"<<std::endl;
			}
	}
}

void compute_EVec(Matrice<double>& H, Write& w){
	Lapack<double> ES(H.ptr(),H.size(),'S');
	Vecteur<double> EVal(H.size());
	ES.eigensystem(EVal);
	w<<H;
}

void compute_EVec(Matrice<std::complex<double> >& H, Write& w){
	Lapack<std::complex<double> > ES(H.ptr(),H.size(),'H');
	Vecteur<double> EVal(H.size());
	ES.eigensystem(EVal);
	w<<H;
}

template<typename M>
void init_nts(Matrice<M> const& U, unsigned int N_n, Write& w){
	unsigned int k(0),N_site(U.size());
	Array2D<unsigned int> nts(N_site*N_n/2,2);
	for(unsigned int i(0); i<N_site;i++){
		for(unsigned int j(i+1); j<N_site;j++){
			if ( std::abs(U(i,j)) > 1e-4){
				nts(k,0) = i;
				nts(k,1) = j;
				k++;
			}
		}
	}
	w<<nts;
}

void check(std::string lattice){
	std::cout<<lattice<<std::endl;
	Read r(lattice.c_str());
	unsigned int N_spin(0), N_m(0), N_n(0);
	bool is_complex;

	r>>is_complex>>N_spin>>N_m>>N_n;
	Array2D<unsigned int> nts(N_spin*N_n*N_m/2,2);
	r>>nts;
	std::cout<<"N_spin="<<N_spin<<" N_m="<<N_m<<" N_n="<<N_n<<std::endl;
	std::cout<<"nts="<<std::endl;
	std::cout<<nts<<std::endl;
	if(is_complex){
		Matrice<std::complex<double> > H;
		r>>H;
		std::cout<<"EVec="<<std::endl;
		std::cout<<H<<std::endl;
	} else {
		Matrice<double> H;
		r>>H;
		std::cout<<"EVec="<<std::endl;
		std::cout<<H<<std::endl;
	}
}
