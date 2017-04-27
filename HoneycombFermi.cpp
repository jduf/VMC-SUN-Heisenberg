#include "HoneycombFermi.hpp"

HoneycombFermi::HoneycombFermi(System const& s):
	System(s),
	Honeycomb<double>(set_ab(),2,"honeycomb-fermi")
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		system_info_.text("HoneycombFermi :");
		system_info_.item("Each color has the same Hamiltonian.");
	}
}

/*{method needed for running*/
void HoneycombFermi::compute_H(){
	H_.set(n_,n_,0);

	double t(-1.0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*t:t);
	}
	H_ += H_.transpose();
}

void HoneycombFermi::create(){
	compute_H();
	diagonalize(true);
	for(unsigned int c(0);c<N_;c++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_(c);j++){
				EVec_[c](i,j) = H_(i,j);
			}
		}
	}
}

Matrix<double> HoneycombFermi::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 1.5;
	tmp(1,0) =-0.5*sqrt(3.0);
	tmp(0,1) = 1.5;
	tmp(1,1) = 0.5*sqrt(3.0);
	return tmp;
}

unsigned int HoneycombFermi::unit_cell_index(Vector<double> const& x) const {
	Vector<double> match(2,0.0);
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 0; }
	match(0) = 1.0/3.0;
	match(1) = 1.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 1; }
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void HoneycombFermi::display_results(){
	compute_H();
	draw_lattice(false,true,false,ref_(3)?(dir_nn_[4]+dir_nn_[3])*1.5:this->dir_nn_[3]*1.25+this->dir_nn_[4]*0.25,"","Fermi");
}

void HoneycombFermi::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="honeycomb-fermi";
	display_results();
}
/*}*/
