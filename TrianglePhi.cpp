#include "TrianglePhi.hpp"

TrianglePhi::TrianglePhi(Container const& param):
	Triangle<std::complex<double> >(param,"triangle-phi")
{}

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

void TrianglePhi::compute_P(Matrix<std::complex<double> >& Px, Matrix<std::complex<double> >& Py){
	Px.set(n_,n_,0.0);
	Py.set(n_,n_,0.0);
	for(unsigned int i(0); i < n_; i++){
		/*horizontal hopping*/
		if( (i % Ly_)  < Ly_ - 1 ){ Px(i,i+1) = 1; }
		else { Px(i,i+1-Lx_) = bc_; }
		/*vertical hopping*/
		if( i+Lx_ < n_ ){ Py(i,i+Lx_) = 1; }
		else { Py(i,i-(Ly_-1)*Lx_) = bc_; }
	}
}

void TrianglePhi::create(double phi){
	phi_=phi;
	compute_T();
	diagonalize_T('H');
	for(unsigned int color(0);color<N_;color++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_;j++){
				EVec_(i+color*n_,j) = T_(i,j);
			}
		}
	}
}

void TrianglePhi::save(){
	filename_ += "-phi" + tostring(phi_);

	Write w(filename_+".jdbin");
	RST rst;
	rst.text("phi-flux : each neighbouring triangle has a flux of opposite sign");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("ref (wave function)",ref_);
	w("n (particles' number)",n_);
	w("N (N of SU(N))",N_);
	w("m (particles per site' number)",m_);
	w("M (particles' number of each color)",M_);
	w("sts (connected sites)",sts_);
	w("phi (flux per triangle)",phi_);
	w("EVec (unitary matrix)",EVec_);
	w("bc (boundary condition)",bc_);
	w("Lx (x-dimension)",Lx_);
	w("Ly (y-dimension)",Ly_);
}

void TrianglePhi::get_param(Container& param){
	GenericSystem<std::complex<double> >::get_param(param);
	param.set("phi_",phi_);
}

void TrianglePhi::study(){
	compute_T();
	//compute_band_structure();
}
