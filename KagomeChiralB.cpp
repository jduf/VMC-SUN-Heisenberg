#include "KagomeChiralB.hpp"

KagomeChiralB::KagomeChiralB(System const& s, double const& phi):
	System(s),
	Kagome<std::complex<double> >(set_ab(),3,"kagome-chiralB"),
	phi_(phi)
{
	if(phi<=4.0){
		if(status_==3){ init_lattice(); }
		if(status_==2){
			init_fermionic();

			system_info_.text("KagomeChiralB :");
			system_info_.item("Each color has the same Hamiltonian.");
			system_info_.item("6 sites per unit cell.");
			system_info_.item("Flux of "+RST::math(my::tostring(phi)+"\\pi/4")+" per plaquette.");
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : the flux per plaquette shouldn't be bigger than pi (phi<=3)"<<std::endl; }
}

/*{method needed for running*/
void KagomeChiralB::compute_H(){
	H_.set(n_,n_,0);

	double t(-1.0);
	double phi(phi_*M_PI/4.0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		switch(obs_[0](i,5)){
			case 0:{ H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*t:t); }break;
			case 1:{ H_(obs_[0](i,0),obs_[0](i,1)) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)==1?phi:0.0)); }break;
			case 2:{ H_(obs_[0](i,0),obs_[0](i,1)) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)==2?phi:0.0)); }break;
		}
	}
	H_ += H_.conjugate_transpose();
}

void KagomeChiralB::create(){
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

void KagomeChiralB::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("phi="+my::tostring(phi_));
		Vector<double> param(1,phi_);

		w.add_to_header()->title(s,'<');
		w<<param;
		w.add_to_header()->add(system_info_.get());
	} else { w<<phi_<<" "; }
}

Matrix<double> KagomeChiralB::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 2.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 1.0;
	tmp(1,1) = sqrt(3.0);
	return tmp;
}

unsigned int KagomeChiralB::unit_cell_index(Vector<double> const& x) const {
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
void KagomeChiralB::display_results(){
	compute_H();
	std::string phi(my::tostring(phi_));
	draw_lattice(true,true,false,(dir_nn_[2]+dir_nn_[3])*0.5,"-d:phi "+phi,RST::math("\\phi")+"="+phi);
}

void KagomeChiralB::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="kagome-chiralB";
	display_results();

	//compute_H();
	//plot_band_structure();
}
/*}*/
