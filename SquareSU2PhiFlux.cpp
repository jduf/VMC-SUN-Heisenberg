#include "SquareSU2PhiFlux.hpp"

SquareSU2PhiFlux::SquareSU2PhiFlux(Parseur& P):
	Square<std::complex<double> >(P,"square-phi-flux"),
	phi_(M_PI*P.get<double>("phi"))
{
	if(!P.status() && N_ == 2){
		if(P.get<bool>("study")){
			compute_T();
			band_structure();
			lattice();
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
				filename_ += "-phi+" + tostring(P.get<double>("phi"));

				save();
			} else {
				std::cerr<<"SquareSU2PhiFlux : degeneate"<<std::endl;
			}
		}
	} else {
		std::cerr<<"SquareSU2PhiFlux : this wavefunction is appropriate for N=2"<<std::endl;
	}
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

void SquareSU2PhiFlux::save(){
	Write w(filename_+".jdbin");
	RST rst;
	rst.text("new test with +- phi flux per plaquette");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("wf (wave function)",wf_);
	w("N (N of SU(N))",N_);
	w("m (m=n/N)",m_);
	w("bc (boundary condition)",bc_);
	w("Lx (x-dimension)",Lx_);
	w("Ly (y-dimension)",Ly_);
	w("phi/pi (phi-flux)",phi_/M_PI);
	w("sts (connected sites)",sts_);
	w("EVec (unitary matrix)",EVec_);
}

void SquareSU2PhiFlux::compute_P(){
	Px_.set(n_,n_,0.0);
	Py_.set(n_,n_,0.0);
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
}

void SquareSU2PhiFlux::band_structure(){
	compute_P();

	//std::cout<<T_*Px_-Px_*T_<<std::endl;
	//std::cout<<T_*Py_-Py_*T_<<std::endl;

	Matrix<std::complex<double> > TP(T_+std::complex<double>(3.)*Px_+std::complex<double>(7.)*Py_);
	Vector<std::complex<double> > eval;
	Matrix<std::complex<double> > evec;
	Lapack<std::complex<double> > ES(&TP,false,'G');
	ES.eigensystem(&eval,&evec);
	Vector<double> kx(n_);
	Vector<double> ky(n_);
	Vector<double> E(n_);
	for(unsigned int i(0);i<n_;i++){
		kx(i) = log(projection(Px_,evec,i,i)).imag()/N_;
		ky(i) = log(projection(Py_,evec,i,i)).imag()-kx(i);
		E(i) = projection(T_,evec,i,i).real();
	}
	save_band_structure(kx,ky,E);
}

void SquareSU2PhiFlux::lattice(){
	PSTricks ps(filename_+"-lattice");
	ps.add("\\begin{pspicture}l(15,15)%"+filename_+"-lattice");
	double prop(0.0);
	std::string color("");
	for(unsigned int i(0);i<sts_.row();i++){
		prop = round(std::abs(T_(sts_(i,0),sts_(i,1)).imag())*7,7);
		if(std::abs(prop)<1e-14){ color = "blue";}
		else{color = "green";}
		switch(H_(sts_(i,0),sts_(i,1))){
			case 1:
				{
					ps.line("->", sts_(i,0)%Lx_, sts_(i,0)/Ly_, sts_(i,1)%Lx_, sts_(i,1)/Ly_, "linewidth="+tostring(prop)+"pt,linecolor="+color);
					break;
				}			
			case -1:
				{
					ps.line("->", sts_(i,0)%Lx_, sts_(i,0)/Ly_,-1, sts_(i,1)/Ly_, "linewidth="+tostring(prop)+"pt,linecolor="+color);
					break;
				}
			case -2:
				{
					ps.line("->", sts_(i,0)%Lx_, sts_(i,0)/Ly_, sts_(i,1)%Lx_, -1, "linewidth="+tostring(prop)+"pt,linecolor="+color);
					break;
				}
			default:
				{
					std::cerr<<"une conditon au bord n'est pas correctement définie"<<std::endl;
				}
		}
	}
	for(unsigned int i(0);i<n_;i++){
		ps.put(i%Lx_,i/Ly_,"\\tiny{"+tostring(i)+"}");
	}
	
	Vector<double> ada(n_,0);
	diagonalize_T('H');

	double r(0.2);
	Vector<double> tmp(2);
	double max(occupation_number(ada));
	for(unsigned int i(0);i<n_;i++){
		tmp(0) = round(ada(i),7);
		tmp(1) = round((max-ada(i))/max,7);
		ps.add("\\rput("+tostring(i%Lx_)+","+tostring(i/Ly_)+"){%");
		ps.pie(tmp,r,"chartColor=color,userColor={blue,white}");
		ps.add("}");
		ps.put(i%Lx_+r*0.5, i/Ly_+r*1.2, "\\tiny{"+tostring(i)+"}");
	}

	ps.frame(-0.5,-0.5,Lx_-0.5,Ly_-0.5,"linecolor=red");
	ps.frame(-0.5,-0.5,1.5,0.5,"linecolor=red,linestyle=dashed");
	ps.add("\\end{pspicture}");
}
