#include "TriangleJastrow.hpp"

TriangleJastrow::TriangleJastrow(Parseur& P):
	Triangle<double>(P,"triangle-Jastrow"),
	nu_(P.get<double>("nu")),
	nn_(n_,z_),
	sl_(n_),
	omega_(N_,N_,1.0)
{
	if(!P.status()){
		if(P.get<bool>("study")){
			lattice();	
		} else {
			filename_ += "-N" + tostring(N_);
			filename_ += "-S" + tostring(n_);
			filename_ += "-" + tostring(Lx_) + "x" + tostring(Ly_);
			if(bc_ == 1){ filename_ += "-P";} 
			else { filename_ += "-A";}
			filename_ += "-nu+" + tostring(nu_);

			compute_nn();
			compute_sublattice();
			compute_omega();
			save();
		}
	} else {
		std::cerr<<"TriangleJastrow : need to provide nu"<<std::endl;
	}
}

TriangleJastrow::~TriangleJastrow(){}

void TriangleJastrow::save(){
	Write w(filename_+".jdbin");
	RST rst;
	rst.text("Staggered magnetic field on the triangle lattice, Becca's idea to mimic an on site chemical potential");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("wf (wave function)",wf_);
	w("N (N of SU(N))",N_);
	w("m (m=n/N)",m_);
	w("bc (boundary condition)",bc_);
	w("Lx (x-dimension)",Lx_);
	w("Ly (y-dimension)",Ly_);
	w("nu (jastrow coefficient)",nu_);
	w("nn (nearst neighbours)",nn_);
	w("sl (sublattice)",sl_);
	w("omega (omega)",omega_);
	w("sts (connected sites)",sts_);
}

void TriangleJastrow::compute_nn(){
	unsigned int k(0);
	for(unsigned int i(0);i<n_;i++){
		k=0;
		for(unsigned int j(0);j<n_;j++){
			if(H_(i,j)){ 
				nn_(i,k) = j;
				k++;
			}
		}
	}
}

void TriangleJastrow::compute_sublattice(){
	unsigned int k(0);
	for(unsigned int i(0);i<n_;i++){
		sl_(i) = k % N_;
		if((i+1)%Lx_==0){
			k++;
		}
		k++;
	}
}

void TriangleJastrow::compute_omega(){
	if(N_==2){
		omega_(1,1) = -1.0;
	}
	if(N_==3){
		omega_(1,1) = std::polar(1.0,2.0*M_PI/3.0);
		omega_(2,2) = std::polar(1.0,2.0*M_PI/3.0);
		omega_(1,2) = std::polar(1.0,4.0*M_PI/3.0);
		omega_(2,1) = std::polar(1.0,4.0*M_PI/3.0);
	}
}

void TriangleJastrow::lattice(){
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

	double r(0.2);
	
	Read reseau("lattice");
	Matrix<unsigned int> ada(n_,N_);
	reseau>>ada;
	double m(ada.max());
	Vector<double> tmp(N_);
	for(unsigned int i(0);i<n_;i++){
		for(unsigned int c(0);c<N_;c++){
			tmp(c) = ada(i,c)/m;
		}
		ps.add("\\rput("+tostring(i%Lx_)+","+tostring(i/Ly_)+"){%");
		ps.pie(tmp,r,"chartColor=color,userColor={blue,red,green}");
		ps.add("}");
		ps.put(i%Lx_+r*0.5, i/Ly_+r*1.2, "\\tiny{"+tostring(i)+"}");
	}

	ps.frame(-0.5,-0.5,Lx_-0.5,Ly_-0.5,"linecolor=red");
	ps.frame(-0.5,-0.5,0.5,0.5,"linecolor=red,linestyle=dashed");
	ps.add("\\end{pspicture}");
}
