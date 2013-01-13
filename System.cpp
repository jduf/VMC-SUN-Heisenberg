#include "System.hpp"

System::System(unsigned int N_spin, unsigned int N_m, unsigned int dim):
	N_spin(N_spin),
	N_m(N_m),
	N_site(N_m*N_spin),
	dim(dim),
	nts(new unsigned int[N_spin*N_m*dim]),
	U(),
	Nx(0),
	Ny(0)
{
	if(dim==1){Nx = N_site; Ny=1;}
	if(dim==2){Nx = 4; Ny = 4;}
	create_U(dim);
	create_nts(dim);
}

System::~System(){
	delete nts;
}

void System::create_U(unsigned int dim){
	arma::Mat<double> T(arma::zeros(N_site,N_site));
	if(dim==1){
		T(0,1)=-1.0;
		T(0,N_site-1)=1.0;
		for(unsigned int i(1); i< N_site-1; i++){
			T(i,i-1) = -1.0;
			T(i,i+1) = -1.0;
		}
		T(N_site-1,0)=1.0;
		T(N_site-1,N_site-2)=-1.0;
	}
	if(dim==2){
		unsigned int j(0);
		for(unsigned int i(0); i< N_site-1; i++){
			if((i+1) % Nx == 0){T(i,j*Nx) = 1;j++;}
			else {T(i,i+1) = -1.0;}
			if(i+Nx < N_site){ T(i,i+Nx) = -1.0;}
			else{ T(i,i+Nx-N_site) = 1.0;}
		}
		T(Nx-1,N_site-1) = 1;
		T(N_site-1,N_site-Ny) = 1;
		T = T+T.t();
	}
	arma::Col<double> EVal;
	arma::eig_sym(EVal,U,T);

	unsigned int val(3);
	arma::Mat<double> A(T);
	//A = A-EVal(val)*arma::eye(N_site,N_site);
	arma::Col<double> X(N_site);
	for(unsigned int i(0);i<N_site;i++){
		X(i) = U(i,val);
	}
	(A*X).print();
}

void System::create_nts(unsigned int dim){
	if(dim==1){
		for(unsigned int i(0);i<N_site;i++){
			nts[i] = (i+1) % N_site;
		}
	} 
	if(dim==2) {
		for(unsigned int j(0);j<Ny;j++){
			for(unsigned int i(0);i<Nx;i++){
				nts[dim*(j*Nx+i)] = (i+1) % Nx + j*Nx;
				nts[dim*(j*Nx+i)+1] = ((j+1)*Nx) % N_site + i;
			}
		}
	}
	//std::cout<<nts[dim*5]<<" "<<nts[dim*5+1]<<std::endl;
}
