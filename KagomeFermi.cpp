#include "KagomeFermi.hpp"

KagomeFermi::KagomeFermi(System const& s):
	System(s),
	Kagome<double>(set_ab(),3,"kagome-fermi")
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		system_info_.text("KagomeFermi :");
		system_info_.item("Each color has the same Hamiltonian.");
	}
}

/*{method needed for running*/
void KagomeFermi::compute_H(){
	H_.set(n_,n_,0);

	double t(-1.0);
	for(unsigned int i(0);i<obs_[0].nlinks(); i++){
		H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*t:t);
	}
	H_ += H_.transpose();
}

void KagomeFermi::create(){
	compute_H();
	diagonalize(true);
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

Matrix<double> KagomeFermi::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 2.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 1.0;
	tmp(1,1) = sqrt(3.0);
	return tmp;
}

unsigned int KagomeFermi::unit_cell_index(Vector<double> const& x) const {
	if(my::are_equal(x(0),0.0,eq_prec_,eq_prec_)){
		if(my::are_equal(x(1),0.0,eq_prec_,eq_prec_)){ return 0; }
		if(my::are_equal(x(1),0.5,eq_prec_,eq_prec_)){ return 1; }
	}
	if(my::are_equal(x(0),0.5,eq_prec_,eq_prec_)){
		if(my::are_equal(x(1),0.0,eq_prec_,eq_prec_)){ return 2; }
	}
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void KagomeFermi::display_results(){
	compute_H();
	draw_lattice(true,true,false,(dir_nn_[4]+dir_nn_[3])*0.5,"","Fermi");
}

void KagomeFermi::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="kagome-fermi";
	display_results();
}
/*}*/
