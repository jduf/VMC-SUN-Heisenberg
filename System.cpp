#include "System.hpp"

/*Constructors and destructor*/
/*{*/
System::System(unsigned int N_spin, unsigned int N_m, unsigned int dim):
	N_spin(N_spin),
	N_m(N_m),
	N_site(N_m*N_spin),
	Nx(0),
	Ny(0),
	dim(dim),
	nts(new unsigned int[N_spin*N_m*dim]),
	U(N_spin*N_m)
{
	if(dim==1){Nx = N_site; Ny=1;}
	if(dim==2){Nx = 4; Ny = 4;}
	create_U(dim);
	create_nts(dim);
}

System::~System(){
	delete[] nts;
}
/*}*/

/*methods that modify the class*/
/*{*/
void System::create_U(unsigned int dim){
	if(dim==1){
		U(0,1)=-1.0;
		U(0,N_site-1)=1.0;
		for(unsigned int i(1); i< N_site-1; i++){
			U(i,i-1) = -1.0;
			U(i,i+1) = -1.0;
		}
		U(N_site-1,0)=1.0;
		U(N_site-1,N_site-2)=-1.0;
	}
	if(dim==2){
		unsigned int j(0);
		for(unsigned int i(0); i< N_site-1; i++){
			if((i+1) % Nx == 0){U(i,j*Nx) = 1;j++;}
			else {U(i,i+1) = -1.0;}
			if(i+Nx < N_site){ U(i,i+Nx) = -1.0;}
			else{ U(i,i+Nx-N_site) = 1.0;}
		}
		U(Nx-1,N_site-1) = 1;
		U(N_site-1,N_site-Ny) = 1;
		//T = T+T.transpose();
		std::cout<<"implémenter transpose"<<std::endl;
	}
	Lapack ES(U.ptr(),U.size(),'S');
	Vecteur EVec(ES.eigensystem());
	EVec.print();
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
/*}*/
