#include "SquarePiFlux.hpp"

SquarePiFlux::SquarePiFlux(System const& s):
	System(s),
	Square<std::complex<double> >(set_ab(),2,"square-piflux")
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		/*This was an attempt to recover results obtained by Wang and
		 * Vishwanath. But as they are using another representation (using
		 * Majorana fermions), the variational wavefunction that they provide
		 * (the one implemented in this class) can't be used to obtain good
		 * result on SU(4) m=1. The energy is way to high.*/
		system_info_.text("SquarePiFlux :");
		system_info_.item("Each color has the same Hamiltonian.");
		system_info_.item("There is a flux of "+RST::math("\\pi") + "per square plaquette.");
	}
}

/*{method needed for running*/
void SquarePiFlux::compute_H(){
	H_.set(n_,n_,0);

	double phi(M_PI/2.0);
	unsigned int s0;
	unsigned int s1;
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		if(obs_[0](i,3)){ H_(s0,s1) = std::polar(obs_[0](i,4)?bc_:1.0,phi)/2.0; }
		else{ H_(s0,s1) = std::polar(obs_[0](i,4)?bc_:1.0,obs_[0](i,5)?-phi:phi)/2.0; }
	}
	H_ += H_.conjugate_transpose();
}

void SquarePiFlux::create(){
	compute_H();
	diagonalize(true);
	status_=1;
	if(status_==1){
		for(unsigned int c(0);c<N_;c++){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = H_(i,j);
				}
			}
		}
	}
}

Matrix<double> SquarePiFlux::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 1;
	tmp(1,0) = 0;
	tmp(0,1) = 0;
	tmp(1,1) = 2;
	return tmp;
}

unsigned int SquarePiFlux::unit_cell_index(Vector<double> const& x) const {
	return my::are_equal(x(1),0.5,eq_prec_,eq_prec_);
}
/*}*/

/*{method needed for checking*/
void SquarePiFlux::display_results(){
	compute_H();
	draw_lattice(true,true,false,(dir_nn_[2]+dir_nn_[3])*0.5,"","pi-flux");
}

void SquarePiFlux::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="square-piflux";
	display_results();

	plot_band_structure();
}
/*}*/
