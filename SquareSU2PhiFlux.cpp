#include "SquareSU2PhiFlux.hpp"

SquareSU2PhiFlux::SquareSU2PhiFlux(Container const& param):
	Square<std::complex<double> >(param,"square-phi-flux"),
	phi_(M_PI*param.get<double>("phi"))
{
	if(bc_ == 1){ filename_ += "-P";} 
	else { filename_ += "-A";}
	filename_ += "-phi+" + tostring(phi_);
	if(N_ == 2){
		compute_T();
		diagonalize_T('H');
		for(unsigned int color(0);color<N_;color++){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_;j++){
					EVec_(i+color*n_,j) = T_(i,j);
				}
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
	w("ref (wave function)",ref_);
	w("n (particles' number)",n_);
	w("N (N of SU(N))",N_);
	w("m (particles per site' number)",m_);
	w("M (particles' number of each color)",M_);
	w("sts (connected sites)",sts_);
	w("phi/pi (phi-flux)",phi_/M_PI);
	w("EVec (unitary matrix)",EVec_);
	w("bc (boundary condition)",bc_);
	w("Lx (x-dimension)",Lx_);
	w("Ly (y-dimension)",Ly_);
}

void SquareSU2PhiFlux::study(){
	compute_T();
	//band_structure();
	//for(unsigned int i(0);i<n_;i++){
	//kx(i) = log(projection(Px_,evec,i,i)).imag()/N_;
	//ky(i) = log(projection(Py_,evec,i,i)).imag()-kx(i);
	//E(i) = projection(T_,evec,i,i).real();
	//}
	lattice();
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

void SquareSU2PhiFlux::lattice(){
	PSTricks ps(filename_+"-lattice");
	ps.add("\\begin{pspicture}l(15,15)%"+filename_+"-lattice");
	double prop(0.0);
	std::string color;
	Vector<unsigned int> neighbourg;
	for(unsigned int i(0);i<n_;i++){
		neighbourg = get_neighbourg(i);
		prop = round(std::abs(T_(i,neighbourg(0)).imag())*7,7);
		if(std::abs(prop)<1e-14){ color = "blue";}
		else{color = "green";}
		if((i+1) % Lx_ ){
			ps.line("->", i%Lx_, i/Ly_, neighbourg(0)%Lx_, neighbourg(0)/Ly_, "linewidth="+tostring(prop)+"pt,linecolor="+color);
		} else {
			ps.line("->", i%Lx_, i/Ly_, i%Lx_+1, neighbourg(0)/Ly_, "linewidth="+tostring(prop)+"pt,linecolor="+color);
		}
		prop = round(std::abs(T_(i,neighbourg(1)).imag())*7,7);
		if(std::abs(prop)<1e-14){ color = "blue";}
		else{color = "green";}
		if( i+Lx_<this->n_){ 
			ps.line("->", i%Lx_, i/Ly_, neighbourg(1)%Lx_, neighbourg(1)/Ly_, "linewidth="+tostring(prop)+"pt,linecolor="+color);
		} else {
			ps.line("->", i%Lx_, i/Ly_, neighbourg(1)%Lx_, i/Ly_+1, "linewidth="+tostring(prop)+"pt,linecolor="+color);
		}
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
