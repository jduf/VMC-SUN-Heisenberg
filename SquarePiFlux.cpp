#include "SquarePiFlux.hpp"

SquarePiFlux::SquarePiFlux(unsigned int N, unsigned int n, unsigned int m, int bc):
	Square<std::complex<double> >(N,n,m,bc,"square-csl")
{
	rst_.text("Chiral spin liquid, with 2pi/N flux per plaquette");
}

SquarePiFlux::~SquarePiFlux(){}

bool SquarePiFlux::create(double x){
	compute_T();
	diagonalize_T('H');
	for(unsigned int spin(0);spin<N_;spin++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_;j++){
				EVec_(i+spin*n_,j) = T_(i,j);
			}
		}
	}
	return !degenerate_;
}

void SquarePiFlux::compute_T(){
	double t(-1.0);
	double phi(2*M_PI/N_);
	for(unsigned int i(0); i< Ly_; i++){
		for(unsigned int j(0); j< Lx_; j++){
			if(j+1 == Lx_){ T_( i*Lx_ , i*Lx_ + j) = bc_*t; }
			else{ T_(i*Lx_ + j , i*Lx_ + j + 1) = t; }
			if(i+1 == Ly_ ){ T_(j, i*Lx_ + j) = std::polar(bc_*t,-((j%N_)+1)*phi); } 
			//if(i+1 == Ly_ ){ T(j, i*Lx_ + j) = std::polar(bc_*t,-(j%N_)*phi); } 
			else{ T_(i*Lx_ + j, (i+1)*Lx_ + j)= std::polar(t,((j%N_)+1)*phi); }
		}
	}
	std::cerr<<"SquarePiFlux : compute_EVec : new use of polar, check that it is correct"<<std::endl;
	std::cerr<<"                            : modified the flux disposition..."<<std::endl;
	T_ += T_.trans_conj(); 
}

void SquarePiFlux::check(){
	compute_T();
	std::cout<<T_.chop(1e-6)<<std::endl;
}

void SquarePiFlux::study(double E, double DeltaE, Vector<double> const& corr, std::string save_in){
	PSTricks ps(save_in,filename_,false);
	ps.add("\\begin{pspicture}(-1,-1)(13,7.75)%"+filename_);
	ps.put(10,7,"$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$ $E\\pm\\Delta E="+tostring(E)+"\\pm"+tostring(DeltaE)+"$");
	std::string color;
	double ll(1.0); //link length
	Matrix<int> nb;
	unsigned int x0;
	unsigned int x1;
	unsigned int y0;
	unsigned int y1;
	for(unsigned int k(0);k<links_.row();k++){
		if(corr(k)>0){color="blue";}
		else{color="red";}
		x0 = (links_(k,0)%Lx_)*ll;
		y0 = (links_(k,0)/Ly_)*ll;
		x1 = (links_(k,1)%Lx_)*ll;
		y1 = (links_(k,1)/Ly_)*ll;
		nb = get_neighbourg(links_(k,0));
		if(nb(0,1) != 1){
			x1 = x0+1;
		}
		if(nb(1,1) != 1){
			y1 = y0+1;
		}
		ps.line("-",x0,y0,x1,y1,"linewidth="+tostring(corr(k))+"pt,linecolor="+color);
	}
	ps.add("\\end{pspicture}");
}

//{//csl for Vishvanath (uses majorana representation)
//for(unsigned int i(0); i< Ly_; i++){
//for(unsigned int j(0); j< Lx_; j++){
//if(j+1 == Lx_){// x hopping
//H(i*Lx_ , i*Lx_ + j) = t;
//if(i % 2 == 0){
//T(i*Lx_ , i*Lx_ + j) = bc_*t;
//} else {
//T(i*Lx_ , i*Lx_ + j) = -bc_*t;
//}
//} else {
//H( i*Lx_ + j , i*Lx_ + j + 1) = t; 
//if(i % 2 == 0){
//T( i*Lx_ + j , i*Lx_ + j + 1) = t; 
//} else {
//T( i*Lx_ + j , i*Lx_ + j + 1) = -t; 
//}
//}
//if(i+1 == Ly_ ){// y hopping
//H(j, i*Lx_ + j) = t;
//T(j, i*Lx_ + j) = bc_*t;
//} else{
//H(i*Lx_ + j, (i+1)*Lx_ + j) = t;
//T(i*Lx_ + j, (i+1)*Lx_ + j) = t;
//}
//}
//}
//H += H.transpose();
//T += T.transpose();
//} 
