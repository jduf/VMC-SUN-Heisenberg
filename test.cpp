#include "Save.hpp"
#include "Chrono.hpp"
#include "Parseur.hpp"
#include "System.hpp"

#include<string>
#include<iostream>

double energie(System& S,unsigned int N_MC);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int N_spin(4);
	unsigned int N_m(6);
	unsigned int N_n(0);
	unsigned int N_MC(1e6);
	std::string filename;
	std::string lattice;

	P.set("N_spin",N_spin);	
	P.set("N_m",N_m);	
	P.set("N_n",N_n);	
	P.set("N_MC",N_MC);	
	P.set("filename",filename);	
	P.set("lattice",lattice);	

	System S(N_spin,N_m,N_n,filename);


	//S.print();
	//S.swap();
	//std::cout<< S.compute_ratio()<<std::endl;
	//S.update_state();
	//S.print();

	//Chrono t;
	//t.tic();
	std::cout<<N_spin<<" "<<N_m<<" "<<N_MC<<" "<<energie(S,N_MC)<<std::endl;
	//t.tac();
	//std::cerr<<t<<" seconde(s)"<<std::endl;
}

double energie(System& S,unsigned int N_MC){
	double ratio(0.0), energie(0.0);
	unsigned int i(0);
	//Save steps("analyse-2d.dat");
	while(i<N_MC){
		S.swap();
		ratio = S.compute_ratio();
		ratio *= ratio;
		if(ratio>1 || (double)rand()/RAND_MAX <ratio){
			i++;
			S.update_state();
			//steps << S << " " << S.det() <<Save::endl;
			for(unsigned int i(0);i<S.N_nts-1;i += 2){
				S.swap(S.nts[i],S.nts[i+1]);
				energie += S.compute_ratio();
			}
		}
	}
	return -energie/(S.N_site * N_MC); // sign(permutation) => ratio det tjrs - ??? 
}


