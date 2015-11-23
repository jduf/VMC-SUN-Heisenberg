#include "TrianglePhi.hpp"

TrianglePhi::TrianglePhi(unsigned int N, unsigned int n, unsigned int m):
	Triangle<std::complex<double> >(N,n,m,"triangle-phi")
{
	rst_.text("phi-flux : each neighbouring triangle has a flux of opposite sign");
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

void TrianglePhi::study(){
	compute_T();
	//compute_band_structure();
}

void TrianglePhi::create(double phi){
	phi_=phi;
	filename_ += "-phi" + tostring(phi_);
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

void TrianglePhi::save(Write& w) const {
	GenericSystem<std::complex<double> >::save(w);
	w("phi (flux per triangle)",phi_);
}

void TrianglePhi::check(){ } 
