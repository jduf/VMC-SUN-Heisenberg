#include "TrianglePhi.hpp"

TrianglePhi::TrianglePhi(Parseur& P):
	Triangle<std::complex<double> >(P,"triangle-phi"),
	phi_(P.get<double>("phi"))
{
	if(!P.status()){
		if(study_system_){
			compute_T();
			compute_band_structure();
		} else {
			compute_T();
			diagonalize_T('H');
			for(unsigned int color(0);color<N_;color++){
				for(unsigned int i(0);i<n_;i++){
					for(unsigned int j(0);j<m_;j++){
						EVec_(i+color*n_,j) = T_(i,j);
					}
				}
			}
			if(successful_){
				filename_ += "-N" + tostring(N_);
				filename_ += "-S" + tostring(n_);
				filename_ += "-" + tostring(Lx_) + "x" + tostring(Ly_);
				if(bc_ == 1){ filename_ += "-P";} 
				else { filename_ += "-A";}
				filename_ += "-phi" + tostring(phi_);
				save();
			} else {
				std::cerr<<"TrianglePhi : degeneate"<<std::endl;
			}
		}
	}
}

TrianglePhi::~TrianglePhi(){}

void TrianglePhi::compute_T(){
	double t(-1.0);
	for(unsigned int i(0); i < n_; i++){
		/*horizontal hopping*/
		if( (i+1) % Lx_ ){ T_(i,i+1) = t;}	
		else { T_(i+1-Lx_,i) = bc_*t;}
		/*vertical hopping*/
		if( i+Lx_ < n_ ){  T_(i,i+Lx_) = t; } 
		else { T_(i-(Ly_-1)*Lx_,i) = bc_*t;}
		/*diagonal hopping*/
		if( (i+1) % Lx_ && i+Lx_ < n_ ){  T_(i,i+Lx_+1) = std::polar(t,phi_); } 
		else {
			if(i+1 < n_ ){
				if( !((i+1) % Lx_) ){ T_(i,i+1) = std::polar(bc_*t,phi_);}/*x jump across boundary*/
				if( i+Lx_ >= n_ ){  T_(i-Lx_*(Ly_-1)+1,i) = std::polar(bc_*t,-phi_); }/*y jump across boundary*/
			} else {
				T_(0,n_-1) = std::polar(bc_*bc_*t,-phi_);
			}
		}
	}
	/*\warning if I take the transpose, the diagonal will be counted twice*/
	T_ += T_.trans_conj();
}

void TrianglePhi::compute_P(){
	Px_.set(n_,n_,0.0);
	Py_.set(n_,n_,0.0);
	for(unsigned int i(0); i < n_; i++){
		/*horizontal hopping*/
		if( (i % Ly_)  < Ly_ - 1 ){ Px_(i,i+1) = 1; }
		else { Px_(i,i+1-Lx_) = bc_; }
		/*vertical hopping*/
		if( i+Lx_ < n_ ){ Py_(i,i+Lx_) = 1; }
		else { Py_(i,i-(Ly_-1)*Lx_) = bc_; }
	}
}

void TrianglePhi::compute_band_structure(){
	compute_P();

	//std::cout<<T_*Px_-Px_*T_<<std::endl;
	//std::cout<<T_*Py_-Py_*T_<<std::endl;

	Matrix<std::complex<double> > TP(T_+std::complex<double>(3.)*Px_+std::complex<double>(7.)*Py_);
	Vector<std::complex<double> > eval;
	Matrix<std::complex<double> > evec;
	Lapack<std::complex<double> > ES(&TP,false,'G');
	ES.eigensystem(&eval,&evec);
	Vector<double> kx(n_,1);
	Vector<double> ky(n_,1);
	Vector<double> E(n_,1);
	for(unsigned int i(0);i<n_;i++){
		kx(i) = log(projection(Px_,evec,i,i)).imag();
		ky(i) = log(projection(Py_,evec,i,i)).imag();
		E(i) = projection(T_,evec,i,i).real();
	}
	save_band_structure(kx,ky,E);
}

void TrianglePhi::save(){
	Write w(filename_+".jdbin");
	RST rst;
	rst.text("phi-flux : each neighbouring triangle has a flux of opposite sign");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("wf (wave function)",wf_);
	w("N (N of SU(N))",N_);
	w("m (number of unit cell)",m_);
	w("bc (boundary condition)",bc_);
	w("Lx (x-dimension)",Lx_);
	w("Ly (y-dimension)",Ly_);
	w("phi (flux per triangle)",phi_);
	w("sts (connected sites)",sts_);
	w("EVec (unitary matrix)",EVec_);
}

