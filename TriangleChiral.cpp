#include "TriangleChiral.hpp"

TriangleChiral::TriangleChiral(System const& s, double const& phi):
	System(s),
	Triangle<std::complex<double> >(set_ab(N_/m_),N_/m_,"triangle-chiral"),
	phi_(phi)
{
	if(1.0*phi_<=N_/m_){
		if(status_==3){ init_lattice(); }
		if(status_==2){
			init_fermionic();

			system_info_.text("TriangleChiral :");
			system_info_.item("Each color has the same Hamiltonian.");
			system_info_.item("There is a flux of "+RST::math(my::tostring(phi)+"\\pi/"+my::tostring(N_/m_)) + "per triangular plaquette.");

			filename_ += "-phi"+my::tostring(phi_);
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : the flux per plaquette shouldn't be bigger than pi"<<std::endl; }
}

/*{method needed for running*/
void TriangleChiral::compute_H(){
	H_.set(n_,n_,0);

	double t(-1.0);
	double phi(phi_*m_*M_PI/N_);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		switch(obs_[0](i,3)){
			case 0:{ H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*t:t); }break;
			case 1:{ H_(obs_[0](i,0),obs_[0](i,1)) = std::polar((obs_[0](i,4)?bc_*t:t),phi*(2.0*obs_[0](i,5)+1)); }break;
			case 2:{ H_(obs_[0](i,0),obs_[0](i,1)) = std::polar((obs_[0](i,4)?bc_*t:t),phi* 2.0*obs_[0](i,5)); }break;
		}
	}
	H_ += H_.conjugate_transpose();
}

void TriangleChiral::create(){
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

void TriangleChiral::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("phi="+my::tostring(phi_));
		Vector<double> param(1,phi_);

		w.add_to_header()->title(s,'<');
		w<<param;
		w.add_to_header()->add(system_info_.get());
	} else { w<<phi_<<" "; }
}

Matrix<double> TriangleChiral::set_ab(unsigned int const& k) const {
	Matrix<double> tmp(2,2);
	if(k == 3 || k == 6){
		tmp(0,0) = k;
		tmp(1,0) = 0.0;
		tmp(0,1) = 0.5;
		tmp(1,1) = sqrt(3.0)/2.0;
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : unknown unit cell"<<std::endl; }
	return tmp;
}

unsigned int TriangleChiral::unit_cell_index(Vector<double> const& x) const {
	switch(spuc_){
		case 3:
			{
				if(my::are_equal(x(0),0.0,    eq_prec_,eq_prec_)){ return 0; }
				if(my::are_equal(x(0),1.0/3.0,eq_prec_,eq_prec_)){ return 1; }
				if(my::are_equal(x(0),2.0/3.0,eq_prec_,eq_prec_)){ return 2; }
			}break;
		case 6:
			{
				if(my::are_equal(x(0),0.0,    eq_prec_,eq_prec_)){ return 0; }
				if(my::are_equal(x(0),1.0/6.0,eq_prec_,eq_prec_)){ return 1; }
				if(my::are_equal(x(0),2.0/6.0,eq_prec_,eq_prec_)){ return 2; }
				if(my::are_equal(x(0),3.0/6.0,eq_prec_,eq_prec_)){ return 3; }
				if(my::are_equal(x(0),4.0/6.0,eq_prec_,eq_prec_)){ return 4; }
				if(my::are_equal(x(0),5.0/6.0,eq_prec_,eq_prec_)){ return 5; }
			}break;
	}
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void TriangleChiral::display_results(){
	compute_H();
	std::string phi(my::tostring(phi_));
	draw_lattice(false,true,false,dir_nn_[3]*1.75+dir_nn_[4]*0.25,"-d:phi "+phi,RST::math("\\phi")+"="+phi);
}

void TriangleChiral::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="triangle-chiral";
	display_results();

	//compute_H();
	//std::cout<<H_<<std::endl;

	//compute_H();
	//plot_band_structure();
}
/*}*/
