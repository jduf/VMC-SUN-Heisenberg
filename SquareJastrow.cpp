#include "SquareJastrow.hpp"

SquareJastrow::SquareJastrow(unsigned int N, unsigned int n, unsigned int m):
	Square<double>(N,n,m,"square-Jastrow"),
	nn_(n_,z_),
	cc_(N_,N_),
	sl_(n_),
	omega_(N_,N_,1.0)
{
	rst_.text("Staggered magnetic field, Becca's idea to mimic an on site chemical potential");
	compute_nn();
	compute_sublattice();
	compute_omega_cc();
}

SquareJastrow::~SquareJastrow(){}

void SquareJastrow::create(double x){
	//nu_.set(z_,x.size()+1);
	//for(unsigned int i(0);i<z_;i++){
		//for(unsigned int j(0);j<Nfreedom_;j++){
			//nu_(i,j) = x(j);
		//}
		//nu_(i,Nfreedom_) = 0;
	//}
}

void SquareJastrow::compute_nn(){
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

void SquareJastrow::compute_sublattice(){
	unsigned int k(0);
	for(unsigned int i(0);i<n_;i++){
		sl_(i) = k % N_;
		if((i+1)%Lx_==0){
			k++;
		}
		k++;
	}
}

void SquareJastrow::compute_omega_cc(){
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

void SquareJastrow::lattice(Matrix<unsigned int> const& lat){
	PSTricks ps(filename_+"-lattice");
	ps.add("\\begin{pspicture}(15,15)%"+filename_+"-lattice");
	Matrix<int> nb;
	double x0, y0, x1, y1;
	std::string color;
	for(unsigned int i(0);i<n_;i++){
		nb = get_neighbourg(i);
		x0 = i%Lx_;
		y0 = i/Ly_;
		y1 = nb(0,0)/Ly_;
		if((i+1) % Lx_ ){
			x1 = nb(0,0)%Lx_;
			color = "black";
		} else {
			x1 = x0 + 1;
			color = "blue";
		}
		ps.line("-", x0, y0, x1, y1 , "linewidth=1pt,linecolor="+color);

		x1 = nb(1,0)%Lx_;
		if( i+Lx_<this->n_){ 
			y1 = nb(1,0)/Ly_;
			color = "black";
		} else {
			y1 = y0 + 1;
			color = "blue";
		}
		ps.line("-", x0, y0, x1, y1, "linewidth=1pt,linecolor="+color);
	}

	double r(0.2);
	Vector<double> pie(N_);
	double m(lat.max());
	for(unsigned int i(0);i<n_;i++){
		for(unsigned int j(0);j<N_;j++){
			pie(j) = lat(i,j)/m;
			if(pie(j) < 1e-4){pie(j) = 0;}
		}
		ps.add("\\rput("+tostring(i%Lx_)+","+tostring(i/Ly_)+"){%");
		ps.pie(pie,r,"chartColor=color,userColor={red,blue}");
		ps.add("}");

		switch(sl_(i)){
			case 0: { color = "red";} break;
			case 1: { color = "blue";} break;
			case 2: { color = "green";} break;
		}
		ps.put(i%Lx_+r*1.2, i/Ly_+r*1.2, "\\tiny{"+tostring(i)+"}");
	}

	ps.frame(-0.5,-0.5,Lx_-0.5,Ly_-0.5,"linecolor=red");
	ps.add("\\end{pspicture}");
}

void SquareJastrow::save(Write& w) const {
	GenericSystem<double>::save(w);
	w("nn (nearst neighbours)",nn_);
	w("cc (to match nu and x)",cc_);
	w("sl (sublattice)",sl_);
	w("omega (omega)",omega_);
}
