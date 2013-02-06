#include "System.hpp"

/*Constructors and destructor*/
/*{*/
System::System(unsigned int N_spin, unsigned int N_m, unsigned int dim):
	N_spin(N_spin),
	N_m(N_m),
	N_site(N_m*N_spin),
	dim(dim),
	nts(new unsigned int[N_spin*N_m*dim]),
	A(new Matrice<double>[N_spin]),
	Ainv(new Matrice<double>[N_spin]),
	s(new unsigned int[N_spin*N_m]),
	wis(new unsigned int[N_spin*N_m]),
	Nx(0),
	Ny(0)
{
	if(dim==1){Nx = N_site; Ny=1;}
	if(dim==2){Nx = 6; Ny = 6;} // Nx & Ny > 2 sinon problème d'antipériodicité
	Matrice<double> U(N_site);
	create_U(U,dim);
	create_nts(dim);
	init_state(U);
	for(unsigned int i(0);i<2;i++){
		w[i] = 0;
		cc[i] = 0;
		mc[i] = 0;
	}
}

System::~System(){
	delete[] nts;
	delete[] s;
	delete[] wis;
	delete[] A;
	delete[] Ainv;
}
/*}*/

/*operators*/
/*{*/
unsigned int System::operator[](unsigned int i) const {
	return wis[i];
}

std::ostream& operator<<(std::ostream& flux, System const& S){
	for(unsigned int i(0);i<S.N_site;i++){
		flux<<S[i]/S.N_m<<" ";
	}
	return flux;
}
/*}*/

/*methods that return something related to the class*/
/*{*/
double System::compute_ratio(){
	if(mc[0] == mc[1]){
		w[0] = -1;
		w[1] = -1;
		return -1;
	} else { 
		w[0] = 0.0;
		w[1] = 0.0;
		for(unsigned int i(0);i<N_m;i++){
			w[0] += Ainv[mc[0]](cc[0],i)*A[mc[1]](i,cc[1]);
			w[1] += Ainv[mc[1]](cc[1],i)*A[mc[0]](i,cc[0]);
		}
		return w[0]*w[1];
	}
}

double System::det(){
	double d(1.0);
	for(unsigned int i(0);i<N_spin;i++){
		Lapack<double> A_(A[i],'G');
		d *= A_.det();
	}
	return d;
}

void System::print(){
	std::cout<<"{";
	for(unsigned int i(0); i<N_spin; i++){
		std::cout<<"{ ";
		for(unsigned int j(0); j<N_m; j++){
			std::cout<<s[i*N_m+j]<<" ";
		}
		std::cout<<"}";
	}
	std::cout<<"}"<<std::endl;
	//for(unsigned int i(0); i<N_site; i++){
		//std::cout<<wis[i]<<" ";
	//}
	//std::cout<<std::endl;
	for(unsigned int i(0);i<2;i++){
		std::cout<<mc[i]<<" "<< cc[i]<<std::endl;
	}
	for(unsigned int i(0);i<N_spin;i++){
		std::cout<<std::endl;
		A[i].print();
	}
	//for(unsigned int i(0);i<N_spin;i++){
		//std::cout<<std::endl;
		//(A[i]*Ainv[i]).print();
	//}
}
/*}*/

/*methods that modify the class*/
/*{*/
void System::create_U(Matrice<double>& U, unsigned int dim){
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
	Vecteur<double> Eval(ES.eigensystem());
	if(fabs(Eval(N_m-1) - Eval(N_m)) < 1e-5){
		std::cerr<<"les valeurs propres sont dégénérées au niveau de Fermi, N_m="<<N_m<<" N_spin="<<N_spin<<std::endl;
		Eval.print();
	}
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
}

void System::init_state(Matrice<double> const& U){
	srand(time(NULL)^(getpid()<<16));
	unsigned int site(0);
	unsigned int N_as(N_site);
	unsigned int available_sites[N_site];
	for(unsigned int i(0); i < N_as; i++){
		available_sites[i]  = i;	
	}

	for(unsigned int i(0); i < N_spin; i++){
		A[i] = Matrice<double> (N_m);
		Ainv[i] = Matrice<double> (N_m);
		for(unsigned int j(0); j < N_m; j++){
			//l = s[i*N_m+j];
			site = rand() % N_as;
			s[j+i*N_m] = available_sites[site];
			wis[available_sites[site]] = i*N_m+j; 	
			for(unsigned int k(0); k < N_m; k++){
				A[i](k,j) = U(available_sites[site],k);
			}
			for(unsigned int k(site); k < N_as-1; k++){
				available_sites[k] = available_sites[k+1];
			}
			N_as--;
		}
		Ainv[i] = A[i];
		Lapack<double> A_(Ainv[i].ptr(),Ainv[i].size(),'G');
		A_.inv();
	}
}

void System::swap() {
	mc[0] = rand() % N_spin;
	mc[1] = rand() % N_spin;
	while(mc[0]==mc[1]){
		mc[1] = rand() % N_spin;
	}
	cc[0] = rand() % N_m;
	cc[1] = rand() % N_m;
}

void System::swap(unsigned int a, unsigned int b) {
	mc[0] = wis[a] / N_m;
	mc[1] = wis[b] / N_m;
	cc[0] = wis[a] % N_m;
	cc[1] = wis[b] % N_m;
}

void System::update_state(){
	unsigned int a(0),b(0);
	a = mc[0]*N_m+cc[0];
	b = mc[1]*N_m+cc[1];

	unsigned int s_tmp(0);
	s_tmp = s[b];
	s[b] = s[a];
	s[a] = s_tmp;

	a = s[a];
	b = s[b];

	s_tmp = wis[b];
	wis[b] = wis[a];
	wis[a] = s_tmp;

	double tmp(0.0);
	for(unsigned int i(0); i<N_m; i++){
		tmp = A[mc[0]](i,cc[0]);
		A[mc[0]](i,cc[0]) = A[mc[1]](i,cc[1]);
		A[mc[1]](i,cc[1]) = tmp;
	}

	//if(mc[0] == mc[1]){
		//Ainv[mc[0]] = A[mc[0]];
		//Lapack A_(Ainv[mc[0]].ptr(),Ainv[mc[0]].size(),'G');
		//A_.inv();
	//} else {
		//Matrice mat(N_m);
		//for(unsigned int m(0);m<2;m++){
			//for(unsigned int i(0);i<N_m;i++){
				//for(unsigned int j(0);j<N_m;j++){
					//tmp = 0.0;
					//for(unsigned int k(0.0);k<N_m;k++){
						//tmp += Ainv[mc[m]](i,k)*A[mc[m]](k,cc[m]);
					//}
					//if(cc[m] == i){
						//tmp -= 1;
					//}
					//tmp /= w[m];
					//mat(i,j) = tmp*Ainv[mc[m]](cc[m],j);
				//}
			//}
			//Ainv[mc[m]] -= mat;
		//}
	//}
	
	if(mc[0] == mc[1]){
		Ainv[mc[0]] = A[mc[0]];
		Lapack<double> A_(Ainv[mc[0]].ptr(),Ainv[mc[0]].size(),'G');
		A_.inv();
	} else {
		for(unsigned int m(0);m<2;m++){
			Ainv[mc[m]] = A[mc[m]];
			Lapack<double> Ainv_(Ainv[mc[m]].ptr(),Ainv[mc[m]].size(),'G');
			Ainv_.inv();
		}
	}
}
/*}*/
