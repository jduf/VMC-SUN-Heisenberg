#include "SquareMu.hpp"

SquareMu::SquareMu(unsigned int N, unsigned int n, unsigned int m):
	Square<double>(N,n,m,"square-mu")
{}

SquareMu::~SquareMu(){}

void SquareMu::compute_T(unsigned int alpha){
	double t(-1.0);
	for(unsigned int i(0); i < n_; i++){
		/*chemical potential*/
		if( (i-alpha) % N_ == 0 && i >= alpha){ T_(i,i) = mu_/2; }
		/*horizontal hopping*/
		if( (i+1) % Lx_ ){ T_(i,i+1) = t;}	
		else { T_(i+1-Lx_,i) = bc_*t; alpha++;}
		/*vertical hopping*/
		if( i+Lx_ < n_ ){  T_(i,i+Lx_) = t; } 
		else { T_(i-(Ly_-1)*Lx_,i) = bc_*t;}
	}
	/*\warning if I take the transpose, the diagonal will be counted twice*/
	T_ += T_.transpose();
}

void SquareMu::compute_P(Matrix<double>& Px, Matrix<double>& Py){
	Px.set(n_,n_,0.0);
	Py.set(n_,n_,0.0);
	for(unsigned int i(0); i < n_; i++){
		/*horizontal hopping*/
		if( (i % Ly_)  < Ly_ - N_ ){Px(i,i+N_) = 1; }
		else{ Px(i,i-Ly_+N_) = bc_; }
		/*vertical hopping*/
		if( i+Lx_ < n_ ){
			if( (i+1) % Lx_ ){Py(i,i+Lx_+1) = 1; }
			else { Py(i,i+1) = bc_;}
		} else {
			if( (i+1) % Lx_ ) {  Py(i,i-(Ly_-1)*Lx_+1) = bc_;}
			else { Py(i,0) = bc_*bc_;}
		}
	}
}

void SquareMu::lattice(){
	PSTricks ps(filename_+"-lattice");
	ps.add("\\begin{pspicture}(15,15)%"+filename_+"-lattice");
	Matrix<int> nb;
	for(unsigned int i(0);i<n_;i++){
		nb = get_neighbourg(i);
		if((i+1) % Lx_ ){
			ps.line("-", i%Lx_, i/Ly_, nb(0,0)%Lx_, nb(0,0)/Ly_, "linewidth=1pt,linecolor=black");
		} else {
			ps.line("-", i%Lx_, i/Ly_, i%Lx_+1, nb(0,0)/Ly_, "linewidth=1pt,linecolor=blue");
		}
		if( i+Lx_<this->n_){ 
			ps.line("-", i%Lx_, i/Ly_, nb(1,0)%Lx_, nb(1,0)/Ly_, "linewidth=1pt,linecolor=black");
		} else {
			ps.line("-", i%Lx_, i/Ly_, nb(1,0)%Lx_, i/Ly_+1, "linewidth=1pt,linecolor=blue");
		}
	}

	diagonalize_T('S');

	double r(0.2);
	Vector<double> ada(n_,0);
	double max(occupation_number(ada));
	Vector<double> tmp(2);
	for(unsigned int i(0);i<n_;i++){
		ps.add("\\rput("+tostring(i%Lx_)+","+tostring(i/Ly_)+"){%");
		ps.pie(tmp,r,"chartColor=color,userColor={blue,white}");
		ps.add("}");
		ps.put(i%Lx_+r*0.5, i/Ly_+r*1.2, "\\tiny{"+tostring(i)+"}");
	}

	ps.frame(-0.5,-0.5,Lx_-0.5,Ly_-0.5,"linecolor=red");
	ps.frame(-0.5,-0.5,N_-0.5,0.5,"linecolor=red,linestyle=dashed");
	ps.add("\\end{pspicture}");
}

void SquareMu::create(double mu){
	mu_ = mu;
	for(unsigned int alpha(0);alpha<N_;alpha++){
		compute_T(alpha);
		diagonalize_T('S');
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_;j++){
				EVec_(i+alpha*n_,j) = T_(i,j);
			}
		}
		T_.set(n_,n_,0.0);
	}
}

void SquareMu::save(Write& w) const {
	GenericSystem<double>::save(w);
	w("mu (chemical potential)",mu_);
}

void SquareMu::check(){
	compute_T(1);
	std::cout<<T_.chop(1e-6)<<std::endl;
}
