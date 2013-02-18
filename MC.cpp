#include "Chrono.hpp"
#include "Read.hpp"
#include "Parseur.hpp"
#include "System.hpp"

#include<string>
#include<iostream>

template<typename T> 
double energie(System<T>& S,unsigned int N_MC);

int main(int argc, char* argv[]){
	std::cerr<<"vérifer abs, fabs, abs(complex)"<<std::endl;
	std::cerr<<"vérifer det et compute_ration pour les complex"<<std::endl;

	Chrono t;
	t.tic();

	double E(0.0);
	Parseur P(argc,argv);
	std::string sysname;
	unsigned int N_spin(0), N_m(0), N_n(0), N_MC(0);
	bool is_complex(false);

	P.set("sysname",sysname);	
	P.set("N_MC",N_MC);	

	Read r(sysname.c_str());
	r>>is_complex>>N_spin>>N_m>>N_n;

	Array2D<unsigned int> sts(N_spin*N_m*N_n/2,2);
	r>>sts;
	if(is_complex){
		Matrice<std::complex<double> > EVec(N_m*N_spin);
		r>>EVec;
		//System<std::complex<double> > S(N_spin,N_m,sts,EVec);
		//E=energie(S,N_MC);
	} else {
		Matrice<double> EVec(N_m*N_spin);
		r>>EVec;
		System<double> S(N_spin,N_m,sts,EVec);
		E=energie(S,N_MC);
	}

	std::cout<<N_spin<<" "<<N_m<<" "<<N_MC<<" "<< E<<std::endl;
	t.tac();
	std::cerr<<t<<" seconde(s)"<<std::endl;
}

template<typename T>
double energie(System<T>& S,unsigned int N_MC){
	double ratio(0.0), energie(0.0);
	unsigned int i(0);
	while(i<N_MC){
		S.swap();
		ratio = S.compute_ratio();
		ratio *= ratio;
		if(ratio>1 || (double)rand()/RAND_MAX <ratio){
			i++;
			S.update_state();
			for(unsigned int i(0);i<S.sts.row();i++){
				S.swap(S.sts(i,0),S.sts(i,1));
				energie += S.compute_ratio();
			}
		}
	}
	return -energie/(S.N_site * N_MC); // sign(permutation) => ratio det tjrs - ??? 
}
