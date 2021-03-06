#include "SquareSU2PhiFlux.hpp"

SquareSU2PhiFlux::SquareSU2PhiFlux(unsigned int N, unsigned int n, unsigned int m):
	Square<std::complex<double> >(N,n,m,"square-phi-flux")
{
	if(N_ != 2){
		std::cerr<<"SquareSU2PhiFlux : this wavefunction is appropriate for N=2"<<std::endl;
	}
	rst_.text("new test with +- phi flux per plaquette");
}

SquareSU2PhiFlux::~SquareSU2PhiFlux(){}

void SquareSU2PhiFlux::compute_T(){
	double t(-1.0);
	bool phase(true);
	for(unsigned int i(0);i<n_;i++){
		/*horizontal hopping*/
		if(phase){
			if( (i+1) % Lx_ ){ T_(i,i+1) = std::polar(t,phi_); phase = false; }	
			else { T_(i+1-Lx_,i) =  std::polar(bc_*t,-phi_); }
		} else {
			if( (i+1) % Lx_ ){ T_(i,i+1) = t; phase = true; }	
			else { T_(i+1-Lx_,i) = bc_*t; }
		}
		/*vertical hopping*/
		if( i+Lx_ < n_ ){  T_(i,i+Lx_) = t; }
		else { T_(i-(Ly_-1)*Lx_,i) = bc_*t;}
	}
	T_ += T_.trans_conj();
}

void SquareSU2PhiFlux::compute_P(Matrix<std::complex<double> >& Px, Matrix<std::complex<double> >& Py){
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

void SquareSU2PhiFlux::save(Write& w) const{

	w("phi/pi (phi-flux)",phi_/M_PI);
}

void SquareSU2PhiFlux::lattice(){
	PSTricks ps(filename_+"-lattice");
	ps.add("\\begin{pspicture}l(15,15)%"+filename_+"-lattice");
	double prop(0.0);
	std::string color;
	Matrix<int> nb;
	for(unsigned int i(0);i<n_;i++){
		nb = get_neighbourg(i);
		if(std::abs(prop)<1e-14){ color = "blue";}
		else{color = "green";}
		if((i+1) % Lx_ ){
			ps.line("->", i%Lx_, i/Ly_, nb(0,0)%Lx_, nb(0,0)/Ly_, "linewidth="+tostring(prop)+"pt,linecolor="+color);
		} else {
			ps.line("->", i%Lx_, i/Ly_, i%Lx_+1, nb(0,0)/Ly_, "linewidth="+tostring(prop)+"pt,linecolor="+color);
		}
		if(std::abs(prop)<1e-14){ color = "blue";}
		else{color = "green";}
		if( i+Lx_<this->n_){
			ps.line("->", i%Lx_, i/Ly_, nb(1,0)%Lx_, nb(1,0)/Ly_, "linewidth="+tostring(prop)+"pt,linecolor="+color);
		} else {
			ps.line("->", i%Lx_, i/Ly_, nb(1,0)%Lx_, i/Ly_+1, "linewidth="+tostring(prop)+"pt,linecolor="+color);
		}
	}

	Vector<double> ada(n_,0);
	diagonalize_T('H');

	double r(0.2);
	Vector<double> tmp(2);
	double max(occupation_number(ada));
	for(unsigned int i(0);i<n_;i++){
		ps.add("\\rput("+tostring(i%Lx_)+","+tostring(i/Ly_)+"){%");
		ps.pie(tmp,r,"chartColor=color,userColor={blue,white}");
		ps.add("}");
		ps.put(i%Lx_+r*0.5, i/Ly_+r*1.2, "\\tiny{"+tostring(i)+"}");
	}

	ps.frame(-0.5,-0.5,Lx_-0.5,Ly_-0.5,"linecolor=red");
	ps.frame(-0.5,-0.5,1.5,0.5,"linecolor=red,linestyle=dashed");
	ps.add("\\end{pspicture}");
}

void SquareSU2PhiFlux::create(double phi){
	filename_ += "-phi+" + tostring(phi_);
	phi_=M_PI*phi;
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

void SquareSU2PhiFlux::check(){
	compute_T();
	std::cout<<T_.chop(1e-6)<<std::endl;
}
