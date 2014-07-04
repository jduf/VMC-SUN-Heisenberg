#include "TriangleFermi.hpp"

TriangleFermi::TriangleFermi(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc):
	System(ref,N,m,n,M,bc),
	Triangle<double>(1,1,1,"triangle-fermi")
{
	if(status_==1){
		init_fermionic();
		compute_T();

		system_info_.text("Fermi : all colors experience the same Hamiltonian");
	}
}

/*{method needed for running*/
void TriangleFermi::compute_T(){
	double t(-1.0);
	Matrix<int> nb;
	for(unsigned int i(0); i < n_; i++){
		nb = get_neighbourg(i);
		for(unsigned int j(0);j<3;j++){
			T_(i,nb(j,0)) = nb(j,1)*t;
		}
	}
	T_ += T_.transpose();
}

void TriangleFermi::create(){
	diagonalize_T();
	for(unsigned int c(0);c<N_;c++){
		if(!is_degenerate(c)){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = T_(i,j);
				}
			}
		}
	}
}
/*}*/

/*{method needed for checking*/
void TriangleFermi::lattice(){
	PSTricks ps("./",filename_+"-lattice");
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
		if((i+1) % Lx_ ){ x1 = nb(0,0)%Lx_; }
		else {
			x1 = x0 + 1;
			color = "blue";
		}
		ps.line("-", x0-y0*e1,y0*e2,x1-y1*e1,y1*e2, "linewidth=1pt,linecolor="+color);

		color = "black";
		if((i+1) % Lx_ ){ x1 = nb(1,0)%Lx_; }
		else {
			x1 = x0 + 1;
			color = "blue";
		}
		if( i+Lx_<this->n_){ y1 = nb(1,0)/Ly_; }
		else {
			y1 = y0 + 1;
			color = "blue";
		}
		ps.line("-", x0-y0*e1,y0*e2,x1-y1*e1,y1*e2, "linewidth=1pt,linecolor="+color);

		color = "black";
		x1=nb(2,0)%Lx_;
		if( i+Lx_<this->n_){ y1 = nb(2,0)/Ly_; }
		else {
			y1 = y0 + 1;
			color = "blue";
		}
		ps.line("-", x0-y0*e1,y0*e2,x1-y1*e1,y1*e2, "linewidth=1pt,linecolor="+color);
	}

	diagonalize_T();

	double r(0.2);
	Vector<double> ada(n_,0);
	//double max(occupation_number(ada));
	Vector<double> tmp(2);
	for(unsigned int i(0);i<n_;i++){
		x0 = i%Lx_;
		y0 = i/Ly_;
		ps.add("\\rput("+tostring(x0-y0*e1)+","+tostring(y0*e2)+"){%");
		ps.pie(tmp,r,"chartColor=color,userColor={blue,white}");
		ps.add("}");
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

void TriangleFermi::check(){
	lattice();
}
/*}*/
