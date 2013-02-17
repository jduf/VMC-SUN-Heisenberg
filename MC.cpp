#include "Chrono.hpp"
#include "Read.hpp"
#include "Parseur.hpp"
#include "System.hpp"

#include<string>
#include<iostream>

double energie(System<double>& S,unsigned int N_MC);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	std::string sysname;
	unsigned int N_spin(0), N_m(0), N_n(0), N_MC(0);

	P.set("sysname",sysname);	
	P.set("N_MC",N_MC);	

	Read r(sysname.c_str());
	r>>N_spin>>N_m>>N_n;
	Matrice<double> U(N_m*N_spin);
	r>>U;
	System<double> S(N_spin,N_m,N_n,U);

	Chrono t;
	t.tic();
	std::cout<<N_spin<<" "<<N_m<<" "<<N_MC<<" "<<energie(S,N_MC)<<std::endl;
	t.tac();
	std::cerr<<t<<" seconde(s)"<<std::endl;
}

double energie(System<double>& S,unsigned int N_MC){
	double ratio(0.0), energie(0.0);
	unsigned int i(0);
	while(i<N_MC){
		S.swap();
		ratio = S.compute_ratio();
		ratio *= ratio;
		if(ratio>1 || (double)rand()/RAND_MAX <ratio){
			i++;
			S.update_state();
			for(unsigned int i(0);i<S.N_nts-1;i += 2){
				S.swap(S.nts[i],S.nts[i+1]);
				energie += S.compute_ratio();
			}
		}
	}
	return -energie/(S.N_site * N_MC); // sign(permutation) => ratio det tjrs - ??? 
}


