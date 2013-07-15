#include"Chain.hpp"

Chain::Chain(Parseur& P):
	CreateSystem<double>(P,2)
{
	if(!P.status()){
		mat_type='S';
		std::string filename("chain-N"+tostring(N_spin) + "-S" + tostring(N_site));
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
		save(filename);
	}
}

Chain::~Chain(){ }

void Chain::compute_T(){
	double t(-1.0);
	H(0, N_site -1 ) = t;
	T(0, N_site -1 ) = bc*t;
	for(unsigned int i(0); i< N_site-1; i++){
		H(i,i+1) = -1;
		T(i,i+1) = t;
	}
	H += H.transpose();
	T += T.transpose();
}

void Chain::save(std::string filename){
	Write w(filename+".jdbin");
	RST rst;
	rst.text("Spin chain, all the hopping parameters are real");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("is_complex",false);
	w("N_spin",N_spin);
	w("N_m",N_m);
	w("sts",sts);
	w("EVec",T);
	w("bc",bc);
}
