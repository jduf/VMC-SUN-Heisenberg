#include "Save.hpp"
#include "State.hpp"
#include "System.hpp"
#include "Chrono.hpp"

double energie(System *S);

int main(){
	unsigned int const N_spin(3);
	unsigned int const N_m(2);
	System oneD(N_spin, N_m, 1);	

	Chrono t;
	t.tic();
	std::cout<<energie(&oneD)<<" ";
	t.tac();
	std::cout<<"en "<<t<<" seconde(s)"<<std::endl;
	//
	//State alpha(&oneD,true);
	//alpha.print();
	//State tmp(&oneD,false);
	//double alpha_det(0.0);
	//double tmp_det(0.0);
//
	//alpha_det=alpha.Det();
//
	//tmp = alpha.swap(1,2);
	//std::cout<<"ratio="<<tmp/alpha<<std::endl;
	//
	//alpha = tmp;
	//alpha.print();
	//tmp_det=alpha.Det();
	//std::cout<<"ratio old="<<tmp_det/alpha_det<<std::endl;
//
//// boucle
	//alpha_det = tmp_det;
//
	//tmp = alpha.swap(1,2);
	//std::cout<<"ratio="<<tmp/alpha<<std::endl;
	//
	//alpha = tmp;
	//alpha.print();
	//tmp_det = alpha.Det();
	//std::cout<<"ratio old="<<tmp_det/alpha_det<<std::endl;
//
//// boucle
	//alpha_det = tmp_det;
//
	//tmp = alpha.swap(1,2);
	//std::cout<<"ratio="<<tmp/alpha<<std::endl;
	//
	//alpha = tmp;
	//alpha.print();
	//tmp_det = alpha.Det();
	//std::cout<<"ratio old="<<tmp_det/alpha_det<<std::endl;


	//double alpha_det(0.0);
	//double tmp_det(0.0);
	//double ratio(0.0);
	//for(unsigned int i(0);i<5;i++){
		//tmp = alpha.swap();
		//ratio = tmp/alpha;
		//alpha_det = alpha.Det();
		//alpha = tmp;
		//tmp_det = alpha.Det();
		//std::cout<<ratio<<" "<<tmp_det/alpha_det<<std::endl;
	//}

	//alpha = tmp;
	//double tmp_det(alpha.Det());
	//std::cout<<tmp_det/alpha_det<<std::endl;
	//std::cout<<"alpha=tmp"<<std::endl;
	//alpha.print();
}

double energie(System *S){
	State alpha(S,true);
	State tmp(S,false);
	double ratio(0.0), energie(0.0);
	unsigned int i(0),NMC(40);
	while(i<NMC){
		tmp = alpha.swap();
		ratio = tmp/alpha;
		ratio *= ratio;
		if(ratio>1 || (double)rand()/RAND_MAX <ratio){
			i++;
			alpha = tmp;
			for(unsigned int j(0);j<S->N_site;j++){
				for(unsigned int d(0);d<S->dim;d++){
					tmp = alpha.swap(j,S->nts[S->dim*j+d]);
					energie += tmp/alpha;
				}
			}
		}
	}
	return -energie/(S->N_site * NMC); // sign(permutation) => ratio det tjrs - ??? 
}


