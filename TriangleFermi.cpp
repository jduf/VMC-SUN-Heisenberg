#include "TriangleFermi.hpp"

TriangleFermi::TriangleFermi(Parseur& P):
	Triangle<double>(P,"triangle-fermi")
{
	if(!P.status()){
		if(P.get<bool>("study")){
			compute_T();
			compute_P();
			band_structure();
			lattice();
		} else {
			compute_T();
			diagonalize_T('S');
			for(unsigned int spin(0);spin<N_;spin++){
				for(unsigned int i(0);i<n_;i++){
					for(unsigned int j(0);j<m_;j++){
						EVec_(i+spin*n_,j) = T_(i,j);
					}
				}
			}
			if(successful_){
				filename_ += "-N" + tostring(N_);
				filename_ += "-S" + tostring(n_);
				filename_ += "-" + tostring(Lx_) + "x" + tostring(Ly_);
				if(bc_ == 1){ filename_ += "-P";} 
				else { filename_ += "-A";}
				save();
			} else {
				std::cerr<<"TriangleFermi : degeneate"<<std::endl;
			}
		}
	}
}

TriangleFermi::~TriangleFermi(){}

void TriangleFermi::compute_T(){
	double t(-1.0);
	for(unsigned int i(0); i < n_; i++){
		/*horizontal hopping*/
		if( (i+1) % Lx_ ){ T_(i,i+1) = t;}	
		else { T_(i+1-Lx_,i) = bc_*t;}
		/*vertical hopping*/
		if( i+Lx_ < n_ ){  T_(i,i+Lx_) = t; } 
		else { T_(i-(Ly_-1)*Lx_,i) = bc_*t;}
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
	T_ += T_.transpose();
}

void TriangleFermi::save(){
	Write w(filename_+".jdbin");
	RST rst;
	rst.text("fermi : all colors experience the same Hamiltonian");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("wf (wave function)",wf_);
	w("N (N of SU(N))",N_);
	w("m (m=n/N)",m_);
	w("bc (boundary condition)",bc_);
	w("Lx (x-dimension)",Lx_);
	w("Ly (y-dimension)",Ly_);
	w("sts (connected sites)",sts_);
	w("EVec (unitary matrix)",EVec_);
}

void TriangleFermi::compute_P(){
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

void TriangleFermi::band_structure(){
	//std::cout<<T_*Px_-Px_*T_<<std::endl;
	//std::cout<<T_*Py_-Py_*T_<<std::endl;

	Matrix<double> TP(T_+3.*Px_+7.*Py_);
	Vector<std::complex<double> > eval;
	Matrix<std::complex<double> > evec;
	Lapack<double> ES(&TP,false,'G');
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

void TriangleFermi::lattice(){
	PSTricks ps(filename_+"-lattice");
	ps.add("\\begin{pspicture}(15,15)%"+filename_+"-lattice");
	std::string color("black");
	double prop(1);
	for(unsigned int i(0);i<sts_.row();i++){
		switch(H_(sts_(i,0),sts_(i,1))){
			case 1:
				{
					ps.line("-", sts_(i,0)%Lx_, sts_(i,0)/Ly_, sts_(i,1)%Lx_, sts_(i,1)/Ly_, "linewidth="+tostring(prop)+"pt,linecolor="+color);
					break;
				}			
			case 2:
				{
					ps.line("-", 0, 0, -1, -1, "linewidth="+tostring(prop)+"pt,linecolor=black");
					break;
				}
			case -1:
				{
					ps.line("-", sts_(i,0)%Lx_, sts_(i,0)/Ly_,-1, sts_(i,1)/Ly_, "linewidth="+tostring(prop)+"pt,linecolor=green");
					break;
				}
			case -2:
				{
					ps.line("-", sts_(i,0)%Lx_, sts_(i,0)/Ly_, sts_(i,1)%Lx_, -1, "linewidth="+tostring(prop)+"pt,linecolor=red");
					break;
				}
			case -3:
				{
					ps.line("-", -1, sts_(i,0)/Ly_, sts_(i,1)%Lx_, sts_(i,1)/Ly_, "linewidth="+tostring(prop)+"pt,linecolor=yellow");
					break;
				}
			case -4:
				{
					ps.line("-", sts_(i,0)%Lx_, sts_(i,0)/Ly_, sts_(i,1)%Lx_, -1, "linewidth="+tostring(prop)+"pt,linecolor=blue");
					break;
				}
			default:
				{
					std::cerr<<"une conditon au bord n'est pas correctement définie"<<std::endl;
				}
		}
	}

	diagonalize_T('S');

	double r(0.2);
	Vector<double> ada(n_,0);
	double max(occupation_number(ada));
	Vector<double> tmp(2);
	for(unsigned int i(0);i<n_;i++){
		tmp(0) = round(ada(i),7);
		tmp(1) = round((max-ada(i))/max,7);
		ps.add("\\rput("+tostring(i%Lx_)+","+tostring(i/Ly_)+"){%");
		ps.pie(tmp,r,"chartColor=color,userColor={blue,white}");
		ps.add("}");
		ps.put(i%Lx_+r*0.5, i/Ly_+r*1.2, "\\tiny{"+tostring(i)+"}");
	}

	ps.frame(-0.5,-0.5,Lx_-0.5,Ly_-0.5,"linecolor=red");
	ps.frame(-0.5,-0.5,0.5,0.5,"linecolor=red,linestyle=dashed");
	ps.add("\\end{pspicture}");
}
