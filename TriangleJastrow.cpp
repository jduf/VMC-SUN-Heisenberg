#include "TriangleJastrow.hpp"

TriangleJastrow::TriangleJastrow(unsigned int N, unsigned int n, unsigned int m):
	Triangle<double>(N,n,m,"triangle-Jastrow"),
	nn_(n_,z_),
	cc_(N_,N_),
	sl_(n_),
	omega_(N_,N_,1.0)
{
	rst_.text("Staggered magnetic field on the triangle lattice, Becca's idea to mimic an on site chemical potential");
	compute_nn();
	compute_sublattice();
	compute_omega_cc();
}

TriangleJastrow::~TriangleJastrow(){}

void TriangleJastrow::create(double x){
	//nu_.set(z_,x.size()+1);
	//for(unsigned int i(0);i<z_;i++){
		//for(unsigned int j(0);j<Nfreedom_;j++){
			//nu_(i,j) = x(j);
		//}
		//nu_(i,Nfreedom_) = 0;
	//}
}

void TriangleJastrow::compute_nn(){
	Matrix<int> nb;
	for(unsigned int i(0);i<n_;i++){
		nb = get_neighbourg(i);
		for(unsigned int j(0);j<z_;j++){
			nn_(i,j) = nb(j,0);
		}
		//unsigned int l(z_);
		//for(unsigned int j(0);j<z_;j++){
		//nb = get_neighbourg(nn_(i,j));
		//for(unsigned int k(j);k<j+2;k++){
		//nn_(i,l) = nb(k%z_,0);
		//l++;
		//}
		//}
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

void TriangleJastrow::compute_omega_cc(){
	if(N_==2){
		omega_(1,1) = -1.0;
		cc_(0,0) = 0;
		cc_(0,1) = 1;
		cc_(1,0) = 1;
		cc_(1,1) = 1;
	}
	if(N_==3){
		omega_(1,1) = std::polar(1.0,2.0*M_PI/3.0);
		omega_(2,2) = std::polar(1.0,2.0*M_PI/3.0);
		omega_(1,2) = std::polar(1.0,4.0*M_PI/3.0);
		omega_(2,1) = std::polar(1.0,4.0*M_PI/3.0);
		cc_(0,0) = 0;
		cc_(0,1) = 1;
		cc_(0,2) = 2;
		cc_(1,0) = 1;
		cc_(1,1) = 3;
		cc_(1,2) = 4;
		cc_(2,0) = 2;
		cc_(2,1) = 4;
		cc_(2,2) = 4;
	}
}

void TriangleJastrow::lattice(Matrix<unsigned int> const& lat){
	PSTricks ps(filename_+"-lattice");
	ps.add("\\begin{pspicture}(15,15)%"+filename_+"-lattice");
	double e1(0.5);
	double e2(sqrt(3.0)/2.0);
	//double e1(0);/*will have exactly the same topology except for the link (0,N)*/
	//double e2(1);
	Matrix<int> nb;
	double x0, y0, x1, y1;
	std::string color;
	for(unsigned int i(0);i<n_;i++){
		nb = get_neighbourg(i);
		x0=i%Lx_;
		y0=i/Ly_;

		color = "black";
		y1 = nb(0,0)/Ly_;
		if((i+1) % Lx_ ){
			x1 = nb(0,0)%Lx_;
		} else {
			x1 = x0 + 1;
			color = "blue";
		}
		ps.line("-", x0-y0*e1,y0*e2,x1-y1*e1,y1*e2, "linewidth=1pt,linecolor="+color);

		color = "black";
		if((i+1) % Lx_ ){
			x1 = nb(1,0)%Lx_;
		} else {
			x1 = x0 + 1;
			color = "blue";
		}
		if( i+Lx_<this->n_){ 
			y1 = nb(1,0)/Ly_;
		} else {
			y1 = y0 + 1;
			color = "blue";
		}
		ps.line("-", x0-y0*e1,y0*e2,x1-y1*e1,y1*e2, "linewidth=1pt,linecolor="+color);

		color = "black";
		x1=nb(2,0)%Lx_;
		if( i+Lx_<this->n_){ 
			y1 = nb(2,0)/Ly_;
		} else {
			y1 = y0 + 1;
			color = "blue";
		}
		ps.line("-", x0-y0*e1,y0*e2,x1-y1*e1,y1*e2, "linewidth=1pt,linecolor="+color);
	}

	double r(0.2);
	Vector<double> pie(N_);
	double m(lat.max());
	for(unsigned int i(0);i<n_;i++){
		x0 = i%Lx_;
		y0 = i/Ly_;
		for(unsigned int j(0);j<N_;j++){
			pie(j) = lat(i,j)/m;
			if(pie(j) < 1e-4){pie(j) = 0;}
		}
		ps.add("\\rput("+tostring(x0-y0*e1)+","+tostring(y0*e2)+"){%");
		ps.pie(pie,r,"chartColor=color,userColor={blue,red,green}");
		ps.add("}");

		switch(sl_(i)){
			case 0: { color = "red";} break;
			case 1: { color = "blue";} break;
			case 2: { color = "green";} break;
		}
		ps.put(x0-y0*e1, y0*e2+r*1.2, "\\color{"+color+"}{\\tiny{"+tostring(i)+"}}");
	}

	Matrix<double> xy(4,2);
	xy(0,0) = -0.5*e1;
	xy(0,1) = -0.5;
	xy(1,0) = Ly_*(-e1)-0.5*e1;
	xy(1,1) = Ly_*e2-0.5;
	xy(2,0) = Lx_-0.5*e1 - Ly_*e1;
	xy(2,1) = Ly_*e2-0.5;;
	xy(3,0) = Lx_-0.5*e1;
	xy(3,1) = -0.5;

	ps.polygon(xy,"linecolor=red");
	ps.add("\\end{pspicture}");
}

void TriangleJastrow::save(Write& w) const {
	GenericSystem<double>::save(w);
	w("nn (nearst neighbours)",nn_);
	w("cc (to match nu and x)",cc_);
	w("sl (sublattice)",sl_);
	w("omega (omega)",omega_);
}

