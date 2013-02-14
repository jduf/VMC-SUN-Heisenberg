#include "System.hpp"

/*Constructors and destructor*/
/*{*/
System::System(unsigned int N_spin, unsigned int N_m, unsigned int N_n, std::string filename):
	N_spin(N_spin),
	N_m(N_m),
	N_n(N_n),
	N_site(N_m*N_spin),
	nts(new unsigned int[2*N_spin*N_m*N_n]),
	A(new Matrice<double>[N_spin]),
	Ainv(new Matrice<double>[N_spin]),
	tmp_mat(N_m,0.0),
	s(new unsigned int[N_spin*N_m]),
	wis(new unsigned int[N_spin*N_m])
{
	srand(time(NULL)^(getpid()<<16));
	Matrice<double> U(N_site);
	Read r(filename, U);
	create_nts(U);
	Lapack<double> ES(U.ptr(),U.size(),'S');
	Vecteur<double> EVal(N_site);
	ES.eigensystem(EVal);
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

/*methods that return something related to the class*/
/*{*/
double System::compute_ratio(){
	if(mc[0] == mc[1]){
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
void System::create_nts(Matrice<double> const& U){
	unsigned int k(0);
	for(unsigned int i(0); i<N_site;i++){
		for(unsigned int j(i+1); j<N_site;j++){
			if ( fabs(U(i,j)) > 1e-4){
				nts[k] = i;
				nts[k+1] = j;
				k+=2;
			}
		}
	}
	for(unsigned int i(0);i<2*N_n*N_site;i += 2){
		std::cout<<nts[i]<<" "<<nts[i+1]<<std::endl;
	}
}

void System::init_state(Matrice<double> const& U){
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
	
	double tmp_start(0.0);
	for(unsigned int m(0);m<2;m++){
		for(unsigned int i(0);i<N_m;i++){
			if(cc[m] == i){ tmp_start = -1.0; }
			else { tmp_start = 0.0;}
			for(unsigned int j(0);j<N_m;j++){
				tmp = tmp_start;
				for(unsigned int k(0);k<N_m;k++){
					tmp += Ainv[mc[m]](i,k)*A[mc[m]](k,cc[m]);
				}
				tmp_mat(i,j) = tmp*Ainv[mc[m]](cc[m],j)/w[m];
			}
		}
		Ainv[mc[m]] -= tmp_mat;
	}
}
/*}*/

std::ostream& operator<<(std::ostream& flux, System const& S){
	for(unsigned int i(0);i<S.N_site;i++){
		flux<<S[i]/S.N_m<<" ";
	}
	return flux;
}
