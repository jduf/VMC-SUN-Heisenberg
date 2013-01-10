#include <armadillo>
#include "State.hpp"
#include "System.hpp"

void find_betas(State const& alpha, State beta[], unsigned int N_site);
double energie(unsigned int N_spin, unsigned int N_m, arma::Mat<double> *U);

int main(){
	unsigned int const N_spin(3);
	unsigned int const N_m(2);
	unsigned int const N_site(N_spin*N_m);

	System oneD(N_spin, N_m, 1);

	//arma::Mat<double> U(arma::zeros(N_site,N_site));
	//U(0,1)=-1.0;
	//U(0,N_site-1)=1.0;
	//for(unsigned int i(1); i< N_site-1; i++){
		//U(i,i-1) = -1.0;
		//U(i,i+1) = -1.0;
	//}
	//U(N_site-1,0)=1.0;
	//U(N_site-1,N_site-2)=-1.0;
	//
	//arma::Col<double> EVal;
	//arma::Mat<double> EVec;
	//arma::eig_sym(EVal,EVec,U);
	//
	//U.print();
	//std::cout<<std::endl;
	//EVec.print();
	//std::cout<<std::endl;
	//EVal.print();
//
	//std::cout<<energie(N_spin,N_m,&EVec)<<std::endl;
}

void find_betas(State const& alpha, State beta[], unsigned int N_permutation){
	for(unsigned int i(0);i<N_permutation;i++){
		beta[i] = alpha.swap(i,(i+1) % N_permutation);
	}
}

double energie(unsigned int N_spin, unsigned int N_m, arma::Mat<double> *U){
	unsigned int const N_permutation(N_spin*N_m); // sera différent en 2D
	State alpha(N_spin,N_m,U);
	State beta[N_permutation];
	State tmp;
	double ratio(0.0), energie(0.0);
	unsigned int i(0),NMC(30);
	while(i<NMC){
		tmp = alpha.swap();
		ratio = (tmp.Det() * tmp.Det()) / (alpha.Det() * alpha.Det());
		if(ratio>1 || rand()/RAND_MAX <ratio){
			alpha = tmp;
			find_betas(alpha,beta,N_permutation); //très différent en 2D
			i++;
			for(unsigned int j(0); j<N_permutation;j++){
				energie += beta[j].Det()/alpha.Det();
			}
		} 
	}
	// sign(permutation) => ratio det tjrs - ??? 
	return -energie/(N_m * N_spin * NMC);
}


