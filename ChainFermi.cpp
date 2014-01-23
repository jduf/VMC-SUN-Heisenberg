#include"ChainFermi.hpp"

ChainFermi::ChainFermi(Parseur& P):
	Chain<double>(P,"chain-fermi")
{
	if(!P.status()){
		filename_ += "-N" + tostring(N_);
		filename_ += "-S" + tostring(n_);
		if(m_ % 2 == 0){ 
			filename_ += "-A";
			bc_ = -1;
		} else {
			filename_ += "-P";
			bc_ = 1;
		}
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
	}
}

ChainFermi::~ChainFermi(){ }

void ChainFermi::compute_T(){
	double t(-1.0);
	T_(0, n_ -1 ) = bc_*t;
	for(unsigned int i(0); i< n_-1; i++){
		T_(i,i+1) = t;
	}
	T_ += T_.transpose();
}

void ChainFermi::compute_P(Matrix<double>& P){
	P.set(n_,n_);
	P(n_ -1,0) = bc_;
	for(unsigned int i(0); i< n_-1; i++){
		P(i,i+1) = 1.0;
	}
}

void ChainFermi::save(){
	Write w(filename_+".jdbin");
	RST rst;
	rst.text("Spin ChainFermi, all the hopping parameters are real");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("ref (wave function)",ref_);
	w("N (N of SU(N))",N_);
	w("m (m=n/N)",m_);
	w("sts (connected sites)",sts_);
	w("EVec (unitary matrix)",EVec_);
	w("bc (boundary condition)",bc_);
}
