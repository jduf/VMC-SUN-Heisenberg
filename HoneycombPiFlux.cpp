#include "HoneycombPiFlux.hpp"

HoneycombPiFlux::HoneycombPiFlux(System const& s):
	System(s),
	Honeycomb<double>(set_ab(),4,"honeycomb-piflux")
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		system_info_.text("HoneycombPiFlux :");
		system_info_.item("Each color has the same Hamiltonian.");
		system_info_.item("pi-flux configuration.");
		system_info_.item("4 sites per unit cell.");
	}
}

/*{method needed for running*/
void HoneycombPiFlux::compute_H(){
	H_.set(n_,n_,0);
	double th(1.0);
	double td(-1.0);

	unsigned int s0(0);
	unsigned int s1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		if(obs_[0](i,5)==2 && obs_[0](i,6)==1 && obs_[0](i,3)==4){ H_(s0,s1) = (obs_[0](i,4)?bc_:1)*td; }
		else { H_(s0,s1) = (obs_[0](i,4)?bc_:1)*th; }
	}
	H_ += H_.transpose();
}

void HoneycombPiFlux::create(){
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

Matrix<double> HoneycombPiFlux::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 3.0;
	tmp(1,0) = 0;
	tmp(0,1) = 0.0;
	tmp(1,1) = sqrt(3.0);
	return tmp;
}

unsigned int HoneycombPiFlux::unit_cell_index(Vector<double> const& x) const {
	Vector<double> match(2,0);
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 0; }
	match(0) = 1.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 1; }
	match(0) = 0.5;
	match(1) = 0.5;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 2; }
	match(0)+= 1.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 3; }
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return 4;
}
/*}*/

/*{method needed for checking*/
void HoneycombPiFlux::display_results(){
	compute_H();
	draw_lattice(true,true,false,ref_(3)?dir_nn_[3]*1.5+(dir_nn_[4]+dir_nn_[5])*0.75:dir_nn_[3]*1.5+(dir_nn_[4]+dir_nn_[5])*0.25,"","Pi-flux");
}

void HoneycombPiFlux::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="honeycomb-piflux";
	display_results();
}
/*}*/
