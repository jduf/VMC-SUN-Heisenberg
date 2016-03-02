#include "TriangleJastrow.hpp"

TriangleJastrow::TriangleJastrow(System const& s, Matrix<double> const& nu):
	System(s),
	Triangle<double>(set_ab(),2,"triangle-jastrow")
{
	if(status_==3){ init_lattice(); }
	init_bosonic(z_,nu);

	system_info_.text("Staggered magnetic field on the triangle lattice, Becca's idea to mimic an on site chemical potential");
	system_info_.item("only works on a two sites per unit cell system (therefore only for SU(2))");

	status_--; 
}

/*{method needed for running*/
void TriangleJastrow::create(){
	unsigned int l(z_/2);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		nn_(obs_[0](i,0),obs_[0](i,3)) = obs_[0](i,1);
		nn_(obs_[0](i,1),l+obs_[0](i,3))=obs_[0](i,0);
		if(!(i%l)){ sl_(obs_[0](i,0)) = obs_[0](i,5); }
	}

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

void TriangleJastrow::save_param(IOFiles& w) const {
	w("nn (nearst neighbours)",nn_);
	w("cc (to match nu and x)",cc_);
	w("sl (sublattice)",sl_);
	w("omega (omega)",omega_);

	GenericSystem<double>::save(w);
}

Matrix<double> TriangleJastrow::set_ab() const {
	Matrix<double> tmp(2,2);
	std::cerr<<"implement a 2 sites unit cell"<<std::endl;
	return tmp;
}

unsigned int TriangleJastrow::unit_cell_index(Vector<double> const& x) const {
	std::cerr<<"implement this on a 2 sites per unit cell"<<std::endl;
	return my::are_equal(x(0),0.5);
}
/*}*/

/*{method needed for checking*/
void TriangleJastrow::check(){
	std::cout<<"void SquareJastrow::check() : nothing to do"<<std::endl;
}
/*}*/
