#include"Chain.hpp"

Chain::Chain(Parseur& P):
	CreateSystem<double>(P,2,"chain")
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
		compute_sts();
		compute_T();
		//compute_band_structure();

		diagonalize_T('S');
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

void Chain::compute_T(){
	double t(-1.0);
	T_(0, n_ -1 ) = bc_*t;
	for(unsigned int i(0); i< n_-1; i++){
		T_(i,i+1) = t;
	}
	T_ += T_.transpose();
}

void Chain::compute_P(Matrix<double>& P){
	P.set(n_,n_);
	P(n_ -1,0) = bc_;
	for(unsigned int i(0); i< n_-1; i++){
		P(i,i+1) = 1.0;
	}
}

void Chain::save(std::string filename){
	Write w(filename+".jdbin");
	RST rst;
	rst.text("Spin chain, all the hopping parameters are real");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("wf (wave function)",wf_);
	w("N (N of SU(N))",N_);
	w("m (m=n/N)",m_);
	w("bc (boundary condition)",bc_);
	w("sts (connected sites)",sts_);
	w("EVec (unitary matrix)",EVec_);
}

Vector<unsigned int> Chain::get_neighbourg(unsigned int i){
	Vector<unsigned int> neighbourg(z_);
	neighbourg(0) = i+1;
	neighbourg(0) = i-1;

	return neighbourg;
}
