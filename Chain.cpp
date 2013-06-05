#include"Chain.hpp"

Chain::Chain(Parseur& P):
	CreateSystem<double>(P,2)
{
	is_complex=false;
	mat_type='S';
	filename = "chain-N"+tostring(N_spin) + "-S" + tostring(N_site);
	if(N_m % 2 == 0){ 
		filename += "-A";
		bc = -1;
	} else {
		filename += "-P";
		bc = 1;
	}
	compute_T();
	compute_sts();
	compute_EVec();
	rst.text("Spin chain, all the hopping parameters are real");
	rst.np();
	save();
}

Chain::~Chain(){ }

void Chain::compute_T(){
	double t(-1.0);
	H(0, N_site -1 ) = t;
	T(0, N_site -1 ) = bc*t;
	for(unsigned int i(0); i< N_site-1; i++){
		H(i,i+1) = t;
		T(i,i+1) = t;
	}
	H += H.transpose();
	T += T.transpose();
}
