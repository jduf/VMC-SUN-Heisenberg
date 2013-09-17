#include"Chain.hpp"

Chain::Chain(Parseur& P):
	CreateSystem<double>(P,2),
	Px(this->N_site,this->N_site)
{
	if(!P.status()){
		std::string filename("chain-N"+tostring(N_spin) + "-S" + tostring(N_site));
		if(N_m % 2 == 0){ 
			filename += "-A";
			bc = -1;
		} else {
			filename += "-P";
			bc = 1;
		}
		compute_H();
		compute_P();
		compute_sts();
		//compute_EVec();
		compute_spectrum();
		for(unsigned int spin(0);spin<N_spin;spin++){
			for(unsigned int i(0);i<N_site;i++){
				for(unsigned int j(0);j<N_m;j++){
					EVec(i+spin*N_site,j) = T(i,j);
				}
			}
		}
		save(filename);
	}
}

Chain::~Chain(){ }

void Chain::compute_H(){
	H(0, N_site -1 ) = 1;
	for(unsigned int i(0); i< N_site-1; i++){
		H(i,i+1) = 1;
	}
	H += H.transpose();
}

void Chain::compute_T(){
	double t(-1.0);
	T(0, N_site -1 ) = bc*t;
	for(unsigned int i(0); i< N_site-1; i++){
		T(i,i+1) = t;
	}
	T += T.transpose();
}

void Chain::compute_P(){
	Px(N_site -1,0) = bc;
	for(unsigned int i(0); i< N_site-1; i++){
		Px(i,i+1) = 1.0;
	}
}

void Chain::compute_EVec(){
	double t(-1.0);
	T(0, N_site -1 ) = bc*t;
	for(unsigned int i(0); i< N_site-1; i++){
		T(i,i+1) = t;
	}
	T += T.transpose();
	diagonalize_EVec('S');
}

void Chain::compute_spectrum(){
	Matrix<double> TP(T+Px);
	Matrix<std::complex<double> > EVal;
	Matrix<std::complex<double> > EVec;
	Lapack<double> ES(&TP,false,'G');
	ES.eigensystem(&EVal,&EVec);
	Matrix<double> k(N_site,1);
	Matrix<double> E(N_site,1);
	for(unsigned int i(0);i<N_site;i++){
		k(i) = log(projection(Px,EVec,i,i)).imag();
		E(i) = projection(T,EVec,i,i).real();
	}
	Gnuplot gp("spectrum","1D");
	gp.save_data("spectrum-tot",k,E);
	gp.save_code();
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
	w("EVec",EVec);
	w("bc",bc);
}
