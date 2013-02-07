#include "System.hpp"

template<>
void System<double>::create_U(Matrice<double>& U, unsigned int dim){
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
		U = U+U.transpose();
		//U.print();
	}
	Lapack<double> ES(U.ptr(),U.size(),'S');
	Vecteur<double> EVal(N_site);
	ES.eigensystem(EVal);
	if(fabs(EVal(N_m-1) - EVal(N_m)) < 1e-5){
		std::cerr<<"les valeurs propres sont dégénérées au niveau de Fermi, N_m="<<N_m<<" N_spin="<<N_spin<<std::endl;
		EVal.print();
	}
}

template<>
void System<std::complex<double> >::create_U(Matrice<std::complex<double> >& U, unsigned int dim){
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
			if((i+1) % Nx == 0){U(j*Nx,i) = 1;j++;} // jump in x direction
			else {U(i,i+1) = -1.0;}
			if(i+Nx < N_site){ U(i,i+Nx) = CC(-1.0, (i % Nx)+1);}
			else{ U(i+Nx-N_site,i) = - CC(1.0, (i % Nx)+1);}
		}
		U(Nx-1,N_site-1) = 1;
		U(N_site-Nx,N_site-1) = 1;
		U = U+U.trans_conj();
		U.print();
	}
	Lapack<std::complex<double> > ES(U.ptr(),U.size(),'S');
	Vecteur<double> EVal(N_site);
	ES.eigensystem(EVal);
	//if(fabs(EVal(N_m-1) - EVal(N_m)) < 1e-5){
		//std::cerr<<"les valeurs propres sont dégénérées au niveau de Fermi, N_m="<<N_m<<" N_spin="<<N_spin<<std::endl;
		////EVal.print();
	//}
}
