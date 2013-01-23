#include "Save.hpp"
#include "State.hpp"
#include "System.hpp"
#include "Chrono.hpp"

double energie(System *S);
double factorial(unsigned int n);

int main(){
	//Chrono t;
	//t.tic();
	unsigned int const N_spin(3);
	unsigned int const N_m(4);
	System oneD(N_spin, N_m, 1);	
	
	State alpha(&oneD);
	State beta(oneD.N_m,oneD.N_spin);

	beta = alpha.swap(1,2);
	alpha.print();
	beta.print();
	
	//std::cout<<energie(&oneD)<<" ";
	//t.tac();
	//std::cout<<"en "<<t<<" seconde(s)"<<std::endl;

	//t.tic();
	//unsigned int const N_spin(4);
	//unsigned int const N_m(4);
	//System two(N_spin, N_m, 2);	
	//std::cout<<energie(&two)<<" ";
	//t.tac();
	//std::cout<<"en "<<t<<" seconde(s)"<<std::endl;
}

double factorial(unsigned int m){
	double n(m);
	return exp(n*log(n)-n+0.5*log(2*M_PI*n)) ;
}

double energie(System *S){
	State alpha(S);
	State beta(S->N_m,S->N_spin);
	State tmp(S->N_m,S->N_spin);
	Save file("det-of-the-chosen-states.dat");

	//std::cout<<"# of states : "<< factorial(S->N_site)/(pow(factorial(S->N_m),S->N_spin))<<std::endl;
	double ratio(0.0), energie(0.0);
	unsigned int i(0),NMC(36);
	while(i<NMC){
		tmp = alpha.swap();
		ratio = (tmp.Det() * tmp.Det()) / (alpha.Det() * alpha.Det());
		if(ratio>1 || (double)rand()/RAND_MAX <ratio){
			i++;
			alpha = tmp;
			file<<alpha.Det()<<Save::endl;
			for(unsigned int j(0);j<S->N_site;j++){
				for(unsigned int d(0);d<S->dim;d++){
					beta = alpha.swap(j,S->nts[S->dim*j+d]);
					energie += beta.Det()/alpha.Det();
				}
			}
		}
	}
	return -energie/(S->N_site * NMC); // sign(permutation) => ratio det tjrs - ??? 
}


