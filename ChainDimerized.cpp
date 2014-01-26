#include"ChainDimerized.hpp"

ChainDimerized::ChainDimerized(Parseur& P):
	Chain<double>(P,"chain-dimerized"),
	delta_(P.get<double>("delta"))
{
	ref_(2) = 1;
	unsigned int pps(2);
	if( !(n_ % 2) ){
		if(!P.status()){
			filename_ += "-delta" + tostring(delta_);
			EVec_.set(N_*n_,pps*m_);
			compute_T();
			//compute_band_structure();

			diagonalize_T('S');
			for(unsigned int spin(0);spin<N_;spin++){
				for(unsigned int i(0);i<n_;i++){
					for(unsigned int j(0);j<pps*m_;j++){
						EVec_(i+spin*n_,j) = T_(i,j);
					}
				}
			}
		}
	} else {
		std::cerr<<"ChainDimerized : need an even number of sites"<<std::endl;
	}
}

ChainDimerized::~ChainDimerized(){ }

void ChainDimerized::compute_T(){
	double t(-1.0);
	for(unsigned int i(0); i< n_-2; i += 2){
		T_(i,i+1) = t+delta_;
		T_(i+1,i+2) = t-delta_;
	}
	T_(n_-2,n_-1) = t+delta_;
	T_(n_-1,0) = bc_*(t-delta_);
	T_ += T_.transpose();
}

void ChainDimerized::compute_P(Matrix<double>& P){
	P.set(n_,n_);
	P(n_ -1,0) = bc_;
	for(unsigned int i(0); i< n_-1; i++){
		P(i,i+1) = 1.0;
	}
}

void ChainDimerized::save(){
	Write w(filename_+".jdbin");
	RST rst;
	rst.text("Spin chain, with different hopping term on the odd and even sites");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("ref (wave function)",ref_);
	w("N (N of SU(N))",N_);
	w("m (m=n/N)",m_);
	w("sts (connected sites)",sts_);
	w("delta (dimerization parameter)",delta_);
	w("EVec (unitary matrix)",EVec_);
	w("bc (boundary condition)",bc_);
}
