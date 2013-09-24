#include"Chain.hpp"

Chain::Chain(Parseur& P):
	CreateSystem<double>(P,2),
	Px_(n_,n_)
{
	if(!P.status()){
		std::string filename("chain-N"+tostring(N_) + "-S" + tostring(n_));
		if(m_ % 2 == 0){ 
			filename += "-A";
			bc_ = -1;
		} else {
			filename += "-P";
			bc_ = 1;
		}
		compute_H();
		compute_P();
		compute_T();
		compute_sts();
		//compute_EVec_();
		compute_band_structure();
		for(unsigned int spin(0);spin<N_;spin++){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<m_;j++){
					EVec_(i+spin*n_,j) = T_(i,j);
				}
			}
		}
		save(filename);
	}
}

Chain::~Chain(){ }

void Chain::compute_H(){
	H_(0, n_ -1 ) = 1;
	for(unsigned int i(0); i< n_-1; i++){
		H_(i,i+1) = 1;
	}
	H_ += H_.transpose();
}

void Chain::compute_T(){
	double t(-1.0);
	T_(0, n_ -1 ) = bc_*t;
	for(unsigned int i(0); i< n_-1; i++){
		T_(i,i+1) = t;
	}
	T_ += T_.transpose();
}

void Chain::compute_P(){
	Px_(n_ -1,0) = bc_;
	for(unsigned int i(0); i< n_-1; i++){
		Px_(i,i+1) = 1.0;
	}
}

void Chain::compute_EVec(){
	double t(-1.0);
	T_(0, n_ -1 ) = bc_*t;
	for(unsigned int i(0); i< n_-1; i++){
		T_(i,i+1) = t;
	}
	T_ += T_.transpose();
	diagonalize_EVec('S');
}

void Chain::compute_band_structure(){
	Matrix<double> TP(T_+Px_);
	Vector<std::complex<double> > eval;
	Matrix<std::complex<double> > evec;
	Lapack<double> ES(&TP,false,'G');
	ES.eigensystem(&eval,&evec);
	Vector<double> k(n_,1);
	Vector<double> E(n_,1);
	for(unsigned int i(0);i<n_;i++){
		k(i) = log(projection(Px_,evec,i,i)).imag();
		E(i) = projection(T_,evec,i,i).real();
	}
	Gnuplot gp("spectrum","1D");
	gp.save_data("spectrum-tot",k,E);
	gp.code(",\\\n");
	Vector<unsigned int> index(E.sort());
	k.sort(index);
	gp.save_data("spectrum-sel",k.range(0,m_),E.range(0,m_));
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
	w("N_",N_);
	w("m_",m_);
	w("sts",sts_);
	w("EVec_",EVec_);
	w("bc_",bc_);
}
