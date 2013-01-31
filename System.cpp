#include "System.hpp"

/*Constructors and destructor*/
/*{*/
System::System(unsigned int N_spin, unsigned int N_m, unsigned int dim):
	N_spin(N_spin),
	N_m(N_m),
	N_site(N_m*N_spin),
	dim(dim),
	nts(new unsigned int[N_spin*N_m*dim]),
	U(N_site),
	A(new Matrice[N_spin]),
	Ainv(new Matrice[N_spin]),
	Nx(0),
	Ny(0)
{
	if(dim==1){Nx = N_site; Ny=1;}
	if(dim==2){Nx = 4; Ny = 4;}
	create_U(dim);
	create_nts(dim);
	for(unsigned int i(0);i<2;i++){
		w[i] = 0;
	}
}

System::~System(){
	delete[] nts;
	delete[] A;
	delete[] Ainv;
}
/*}*/

/*methods that return something related to the class*/
/*{*/
void System::update_matrices(unsigned int const mc[], unsigned int const cc[]){
	double tmp(0.0);
	for(unsigned int i(0); i<N_m; i++){
		tmp = A[mc[0]](i,cc[0]);
		A[mc[0]](i,cc[0]) = A[mc[1]](i,cc[1]);
		A[mc[1]](i,cc[1]) = tmp;
	}

	std::cout<<w[0]<<" "<<w[1]<<std::endl;
	Vecteur v1(N_m);
	Vecteur v2(N_m);
	for(unsigned int m(0);m<2;m++){
		for(unsigned int i(0);i<N_m;i++){
			v1(i) = 0.0;
			v2(i) = Ainv[mc[m]](cc[m],i);
			for(unsigned int j(0);j<N_m;j++){
				v1(i) += Ainv[mc[m]](i,j)*A[mc[m]](j,cc[m])/w[m];
			}
			v1(cc[m]) -= 1.0/w[m];
		}
		Ainv[mc[m]] -= v1^v2;
	}

	
	//for(unsigned int m(0);m<2;m++){
		//Ainv[mc[m]] = A[mc[m]];
		//Lapack A_(Ainv[mc[m]].ptr(),Ainv[mc[m]].size(),'G');
		//A_.inv();
		//Ainv[mc[m]].print();
		//std::cout<<std::endl;
	//}
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
		U = U+U.transpose();
	}
	Lapack ES(U.ptr(),U.size(),'S');
	ES.eigensystem();
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
