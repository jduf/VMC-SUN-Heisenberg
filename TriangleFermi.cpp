#include "TriangleFermi.hpp"

TriangleFermi::TriangleFermi(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc):
	System(ref,N,m,n,M,bc),
	Triangle<double>(set_ab(),1,"triangle-fermi")
{
	if(status_==2){
		init_fermionic();

		system_info_.text("Fermi : all colors experience the same Hamiltonian");
		check_lattice();
	}
}

/*{method needed for running*/
void TriangleFermi::compute_H(){
	double t(-1.0);
	Matrix<int> nb;
	for(unsigned int i(0); i < n_; i++){
		nb = get_neighbourg(i);
		for(unsigned int j(0);j<3;j++){
			H_(i,nb(j,0)) = nb(j,1)*t;
		}
	}
	H_ += H_.transpose();
}

void TriangleFermi::create(){
	compute_H();
	diagonalize(false);
	for(unsigned int c(0);c<N_;c++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_(c);j++){
				EVec_[c](i,j) = H_(i,j);
			}
		}
	}
}

unsigned int TriangleFermi::match_pos_in_ab(Vector<double> const& x) const { 
	(void)(x); 
	return 0;
}

Matrix<double> TriangleFermi::set_ab(){
	Matrix<double> tmp(2,2);
	tmp(0,0) = 1;
	tmp(1,0) = 0;
	tmp(0,1) = 0;
	tmp(1,1) = 1;
	return tmp;
}
/*}*/

/*{method needed for checking*/
void TriangleFermi::lattice(){
	Matrix<double> e(2,2);
	e(0,0) = 1.0/3.0;
	e(1,0) = 1.0/3.0;
	e(0,1) = -1.0/2.0;
	e(1,1) = 1.0/2.0;
	Matrix<double> inv_e(2,2);
	inv_e(0,0) = e(1,1);
	inv_e(1,0) =-e(1,0);
	inv_e(0,1) =-e(0,1);
	inv_e(1,1) = e(0,0);
	inv_e/=(e(0,0)*e(1,1)-e(1,0)*e(0,1));

	Matrix<int> nb;
	std::string color("black");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps("./","lattice");
	ps.add("\\begin{pspicture}(-9,-10)(16,10)%"+filename_);
	for(unsigned int i(0);i<n_;i++) {
		xy0 = get_pos_in_lattice(i);
		set_pos_LxLy(xy0);
		set_in_basis(xy0);
		xy0 = (inv_e*LxLy_*xy0).chop();
		nb = get_neighbourg(i);
		ps.put(xy0(0)-0.20,xy0(1)+0.15,my::tostring(i));

		if(nb(0,1)<0){
			color = "red";
			xy1 = xy0;
			xy1(0) += 1.5;
			xy1(1) -= 1.0;
			ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(nb(0,0)));
		} else {
			color = "black";
			if(i+1 == xloop_){
				xy1 = xy0;
				xy1(0) += 1.5;
				xy1(1) -= 1.0;
			} else {
				xy1 = get_pos_in_lattice(nb(0,0));
				set_pos_LxLy(xy1);
				set_in_basis(xy1);
				xy1 = inv_e*LxLy_*xy1;
			}
		}
		xy1 = xy1.chop();
		/*x-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth=1pt,linecolor="+color);

		if(nb(1,1)<0){
			color = "red";
			xy1 = xy0;
			xy1(0) += 1.5;
			xy1(1) += 1.0;
			ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(nb(1,0)));
		} else {
			color = "black";
			xy1 = get_pos_in_lattice(nb(1,0));
			set_pos_LxLy(xy1);
			set_in_basis(xy1);
			xy1 = inv_e*LxLy_*xy1;
		}
		xy1 = xy1.chop();
		/*y-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth=1pt,linecolor="+color);

		if(nb(2,1)<0){
			color = "red";
			xy1 = xy0;
			xy1(0) += 0.0;
			xy1(1) += 2.0;
			ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(nb(2,0)));
		} else {
			color = "black";
			xy1 = get_pos_in_lattice(nb(2,0));
			set_pos_LxLy(xy1);
			set_in_basis(xy1);
			xy1 = inv_e*LxLy_*xy1;
		}
		xy1 = xy1.chop();
		/*xy-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth=1pt,linecolor="+color);
	}

	Vector<double> Lx(2);
	Lx(0) = LxLy_(0,0);
	Lx(1) = LxLy_(1,0);
	Lx = inv_e*Lx;
	Vector<double> Ly(2);
	Ly(0) = LxLy_(0,1);
	Ly(1) = LxLy_(1,1);
	Ly = inv_e*Ly;

	Matrix<double> polygon(4,2);
	polygon(0,0)=0;
	polygon(0,1)=0;
	polygon(1,0)=Lx(0);
	polygon(1,1)=Lx(1);
	polygon(2,0)=Lx(0)+Ly(0);
	polygon(2,1)=Lx(1)+Ly(1);
	polygon(3,0)=Ly(0);
	polygon(3,1)=Ly(1);
	for(unsigned int i(0);i<polygon.row();i++){
		polygon(i,0) -= 0.3; 
		polygon(i,1) += 0.2; 
	}
	ps.polygon(polygon,"linecolor=green");

	ps.add("\\end{pspicture}");
	ps.save(true,true,true);
}

void TriangleFermi::check(){
	//Matrix<int> nb;
	//std::cout<<"######################"<<std::endl;
	//for(unsigned int i(0);i<n_;i++){
	//nb = get_neighbourg(i);
	//std::cout<<"i="<<i<<std::endl;
	//std::cout<<nb<<std::endl;
	//}
	//std::cout<<"######################"<<std::endl;
	//nb = get_neighbourg(3);
	//std::cout<<nb<<std::endl;
	lattice();
}
/*}*/
