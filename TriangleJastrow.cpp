#include "TriangleJastrow.hpp"

TriangleJastrow::TriangleJastrow(Parseur& P):
	Triangle<double>(P,"triangle-Jastrow"),
	nn_(n_,3*z_),
	sl_(n_),
	omega_(N_,N_,1.0)
{
	if(!P.status()){
		compute_nn();
		compute_sublattice();
		compute_omega();
		if(P.get<bool>("study")){
			lattice();	
			std::cout<<nn_<<std::endl;
		} else {
			filename_ += "-N" + tostring(N_);
			filename_ += "-S" + tostring(n_);
			filename_ += "-" + tostring(Lx_) + "x" + tostring(Ly_);
			if(bc_ == 1){ filename_ += "-P";} 
			else { filename_ += "-A";}
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
	w("nn (nearst neighbours)",nn_);
	w("sl (sublattice)",sl_);
	w("omega (omega)",omega_);
	w("sts (connected sites)",sts_);
}

Vector<unsigned int> TriangleJastrow::get_neighbourg(unsigned int i){
	Vector<unsigned int> neighbourg(z_);
	/*0 neighbour*/
	if((i+1)%Lx_){ neighbourg(0) = i+1; } 
	else { neighbourg(0) = (i/Lx_)*Lx_; }
	/*pi/3 neighbour*/
	if((i+1)%Lx_ && i+Lx_<n_){ neighbourg(1) = i+Lx_+1; } 
	else {
		if(i+1<n_){
			if((i+1)%Lx_){ neighbourg(1) = i-Lx_*(Ly_-1)+1; }
			if(i+Lx_<n_ ){ neighbourg(1) = i+1; }
		} else { neighbourg(1) = 0 ;}
	}
	/*2pi/3 neighbour*/
	if(i+Lx_<n_){ neighbourg(2) = i+Lx_; }
	else { neighbourg(2) = i-n_+Lx_; }
	/*pi neighbour*/
	if(i%Lx_){ neighbourg(3) = i-1; }
	else { neighbourg(3) = i+Lx_-1; }
	/*4pi/3 neighbour*/
	if(i%Lx_ && i>=Lx_){ neighbourg(4) = i-Lx_-1; } 
	else {
		if(i!=0){
			if(i<Lx_){ neighbourg(4) = n_-Lx_+i-1;}
			else{ neighbourg(4) = i-1; }
		} else { neighbourg(4) = n_-1 ;}
	}
	/*5pi/3 neighbour*/
	if(i>=Lx_){ neighbourg(5) = i-Lx_; }
	else { neighbourg(5) = n_-Lx_+i; }

	return neighbourg;
}

void TriangleJastrow::compute_nn(){
	Vector<unsigned int> neighbourg;
	for(unsigned int i(0);i<n_;i++){
		neighbourg = get_neighbourg(i);
		for(unsigned int j(0);j<z_;j++){
			nn_(i,j) = neighbourg(j);
		}
		unsigned int l(z_);
		for(unsigned int j(0);j<z_;j++){
			neighbourg = get_neighbourg(nn_(i,j));
			for(unsigned int k(j);k<j+2;k++){
				nn_(i,l) = neighbourg(k%z_);
				l++;
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

void TriangleJastrow::properties(Container& c){
	c.set("N",N_);
	c.set("m",m_);
	c.set("bc",bc_);
	c.set("Lx",Lx_);
	c.set("Ly",Ly_);
	c.set("sts",sts_);
	c.set("nn",nn_);
	c.set("sl",sl_);
	c.set("omega",omega_);
}

void TriangleJastrow::lattice(){
	PSTricks ps(filename_+"-lattice");
	ps.add("\\begin{pspicture}(15,15)%"+filename_+"-lattice");
	double prop(1);
	double e1(0.5);
	double e2(sqrt(3.0)/2.0);
	//double e1(0);/*will have exactly the same topology except for the link (0,N)*/
	//double e2(1);
	for(unsigned int i(0);i<sts_.row();i++){
		double x0(sts_(i,0)%Lx_);
		double y0(sts_(i,0)/Ly_);
		double x1(sts_(i,1)%Lx_);
		double y1(sts_(i,1)/Ly_);
		switch(H_(sts_(i,0),sts_(i,1))){
			case 1:
				{
					ps.line("-", x0-y0*e1,y0*e2,x1-y1*e1,y1*e2, "linewidth="+tostring(prop)+"pt,linecolor=black");
				}break;
			case 2:
				{
					ps.line("-", 0, 0, -e1, -e2, "linewidth="+tostring(prop)+"pt,linecolor=black");
				}break;
			case -1:
				{
					ps.line("-", x0-y0*e1,y0*e2,-1-y1*e1,y1*e2, "linewidth="+tostring(prop)+"pt,linecolor=green");
				}break;
			case -2:
				{
					ps.line("-", x0,0,x1+e1,-e2, "linewidth="+tostring(prop)+"pt,linecolor=red");
				}break;
			case -3:
				{
					ps.line("-", -1-y0*e1,y0*e2,x1-y1*e1,y1*e2, "linewidth="+tostring(prop)+"pt,linecolor=yellow");
				}break;
			case -4:
				{
					ps.line("-", x0,0,x1+e1,-e2, "linewidth="+tostring(prop)+"pt,linecolor=blue");
				}break;
			default:
				{
					std::cerr<<"une conditon au bord n'est pas correctement définie"<<std::endl;
				}
		}
	}

	double r(0.2);
	//Read reseau("lattice");
	//Matrix<unsigned int> ada(n_,N_);
	//reseau>>ada;
	//double m(ada.max());
	//Vector<double> tmp(N_);
	std::string color;
	for(unsigned int i(0);i<n_;i++){
		double x0(i%Lx_);
		double y0(i/Ly_);
		//for(unsigned int c(0);c<N_;c++){
		//tmp(c) = ada(i,c)/m;
		//}
		//ps.add("\\rput("+tostring(i%Lx_)+","+tostring(i/Ly_)+"){%");
		//ps.pie(tmp,r,"chartColor=color,userColor={blue,red,green}");
		//ps.add("}");

		switch(sl_(i)){
			case 0: { color = "red";} break;
			case 1: { color = "blue";} break;
			case 2: { color = "green";} break;
		}

		ps.put(x0-y0*e1, y0*e2+r*1.2, "\\color{"+color+"}{\\tiny{"+tostring(i)+"}}");
	}

	std::cout<<sl_<<std::endl;

	//ps.frame(-0.5,-0.5,Lx_-0.5,Ly_-0.5,"linecolor=red");
	//ps.frame(-0.5,-0.5,0.5,0.5,"linecolor=red,linestyle=dashed");
	ps.add("\\end{pspicture}");
}
