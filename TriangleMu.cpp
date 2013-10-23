#include "TriangleMu.hpp"

TriangleMu::TriangleMu(Parseur& P):
	Triangle<double>(P),
	mu_(P.get<double>("mu"))
{
	if(!P.status()){
		for(unsigned int color(0);color<N_;color++){
			compute_T(color);
			diagonalize_T('S');
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<m_;j++){
					EVec_(i+color*n_,j) = T_(i,j);
				}
			}
			T_.set(n_,n_,0.0);
		}
		if(successful_){
			std::string filename("triangle-stripe");
			filename += "-N" + tostring(N_);
			filename += "-S" + tostring(n_);
			filename += "-" + tostring(Lx_) + "x" + tostring(Ly_);
			if(bc_ == 1){ filename += "-P";} 
			else { filename += "-A";}
			filename += "-mu" + tostring(mu_);
			save(filename);
		} else {
			std::cerr<<"TriangleMu : degeneate"<<std::endl;
		}
	}
}

TriangleMu::~TriangleMu(){}

void TriangleMu::compute_T(unsigned int color){
	double t(-1.0);
	for(unsigned int i(0); i < n_; i++){
		/*chemical potential*/
		if( (i-color) % N_ == 0 && i >= color){ T_(i,i) = mu_/2; }
		/*horizontal hopping*/
		if( (i+1) % Lx_ ){ T_(i,i+1) = t;}	
		else { T_(i+1-Lx_,i) = bc_*t;}
		/*vertical hopping*/
		if( i+Lx_ < n_ ){  T_(i,i+Lx_) = t; } 
		else { T_(i-(Ly_-1)*Lx_,i) = bc_*t; color++;}
		/*diagonal hopping*/
		if( (i+1) % Lx_ && i+Lx_ < n_ ){  T_(i,i+Lx_+1) = t; } 
		else {
			if(i+1 < n_ ){
				if( !((i+1) % Lx_) ){ T_(i,i+1) = bc_*t;}/*x jump across boundary*/
				if( i+Lx_ >= n_ ){  T_(i-Lx_*(Ly_-1)+1,i) = bc_*t; }/*y jump across boundary*/
			} else {
				T_(0,n_-1) = bc_*bc_*t;
			}
		}
	}
	/*\warning if I take the transpose, the diagonal will be counted twice*/
	T_ += T_.transpose();
}

void TriangleMu::save(std::string filename){
	Write w(filename+".jdbin");
	RST rst;
	rst.text("Stripe order : each color lives on its own sublattice");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("wf (wave function)",wf_);
	w("N (N of SU(N))",N_);
	w("m (number of unit cell)",m_);
	w("bc (boundary condition)",bc_);
	w("Lx (x-dimension)",Lx_);
	w("Ly (y-dimension)",Ly_);
	w("mu (chemical potential)",mu_);
	w("sts (connected sites)",sts_);
	w("EVec (unitary matrix)",EVec_);
}

