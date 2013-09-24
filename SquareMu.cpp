#include "SquareMu.hpp"

SquareMu::SquareMu(Parseur& P):
	Square<double>(P),
	mu_(P.get<double>("mu"))
{
	if(!P.status()){
		if(n_==Ly_*Lx_){
			compute_T(0);
			compute_P();	
			compute_band_structure();	
			//for(unsigned int spin(0);spin<N_;spin++){
				//compute_EVec(spin);
				//for(unsigned int i(0);i<n_;i++){
					//for(unsigned int j(0);j<m_;j++){
						//EVec_(i+spin*n_,j) = T_(i,j);
					//}
				//}
			//}
			if(successful_){
				std::string filename("square-stripe");
				filename += "-N" + tostring(N_);
				filename += "-S" + tostring(n_);
				filename += "-" + tostring(Ly_) + "x" + tostring(Lx_);
				if(bc_ == 1){ filename += "-P";} 
				else { filename += "-A";}
				filename += "-mu_" + tostring(mu_);
				save(filename);
			} else {
				std::cerr<<"SquareMu : degeneate"<<std::endl;
			}
		} else {
			std::cerr<<"SquareMu : the cluster is not a square"<<std::endl;
		}
	}
}

SquareMu::~SquareMu(){}

void SquareMu::compute_T(unsigned int spin){
	T_.set(n_,n_,0.0);
	double t(-1.0);
	for(unsigned int i(0); i < n_; i++){
		/*chemical potential*/
		if( (i-spin) % N_ == 0 && i >= spin){ T_(i,i) = mu_/2; }
		/*horizontal hopping*/
		if( (i+1) % Lx_ ){ T_(i,i+1) = t;}	
		else { T_(i+1-Lx_,i) = bc_*t;  spin++; }
		/*vertical hopping*/
		if( i+Lx_ < n_ ){  T_(i,i+Lx_) = t; } 
		else { T_(i-(Ly_-1)*Lx_,i) = bc_*t;}
	}
	/*\warning if I take the transpose, the diagonal will be counted twice*/
	T_ += T_.transpose();
	//show(T,spin%N_+1);
	//diagonalize_EVec('S');
}

void SquareMu::compute_P(){
	for(unsigned int i(0); i < n_; i++){
		/*horizontal hopping*/
		if( (i % Ly_)  < Ly_ - N_ ){Px_(i,i+N_) = 1; }
		else{ Px_(i,i-Ly_+N_) = bc_; }
		/*vertical hopping*/
		if( i+Lx_ < n_ ){
			if( (i+1) % Lx_ ){Py_(i,i+Lx_+1) = 1; }
			else { Py_(i,i+1) = bc_;}
		} else {
			if( (i+1) % Lx_ ) {  Py_(i,i-(Ly_-1)*Lx_+1) = bc_;}
			else { Py_(i,0) = bc_*bc_;}
		}
	}
	std::cout<<T_*Px_-Px_*T_<<std::endl;
	std::cout<<T_*Py_-Py_*T_<<std::endl;
}

void SquareMu::save(std::string filename){
	Write w(filename+".jdbin");
	RST rst;
	rst.text("Stripe order : each color lives on its own sublattice");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("is_complex",false);
	w("N_",N_);
	w("m_",m_);
	w("sts",sts_);
	w("EVec",EVec_);
	w("bc_",bc_);
	w("Ly_",Ly_);
	w("Lx_",Lx_);
	w("mu_",mu_);
}

void SquareMu::show(Matrix<double> const& T,unsigned int spin){
	std::cout<<"T="<<std::endl;
	std::cout<<T<<std::endl;
	std::cout<<"favored sites :"<<std::endl;
	for(unsigned int i(0);i<Ly_;i++){
		for(unsigned int j(0);j<Lx_;j++){
			if(T(i+j*Ly_,i+j*Ly_)!=0){
				std::cout<<spin<<" ";
			} else { 
				std::cout<<0<<" ";
			}
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}
