#include "Squarek2MuPi.hpp"

Squarek2MuPi::Squarek2MuPi(System const& s, double const& mu):
	System(s),
	Square<double>(set_ab(),N_/m_,"square-k2mupi"),
	mu_(mu)
{
	if(spuc_ == 2 && status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();
		same_wf_ = false;

		system_info_.text("Squarek2MuPi :");
		system_info_.item("Néel-like state.");
		system_info_.item("Each color has a different Hamiltonian.");
		system_info_.item(RST::math("\\pi")+"-flux in each plaquette.");

		filename_ += "-mu"+std::string(mu_>=0?"+":"")+my::tostring(mu_);
	}
}

/*{method needed for running*/
void Squarek2MuPi::compute_H(unsigned int const& c){
	H_.set(n_,n_,0);

	double t(-1.0);
	unsigned int s0(0);
	unsigned int s1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		H_(s0,s1) = (obs_[0](i,5) && obs_[0](i,3)?-1.0:1.0)*(obs_[0](i,4)?bc_*t:t);
		if((unsigned int)(obs_[0](i,5))==c%spuc_){ H_(s0,s0) = mu_/2; }
	}
	H_ += H_.transpose();
}

void Squarek2MuPi::create(){
	for(unsigned int c(0);c<N_;c++){
		status_ = 2;
		compute_H(c);
		diagonalize(true);
		if(status_==1){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = H_(i,j);
				}
			}
		} else { c = N_; }
	}
}

void Squarek2MuPi::save_param(IOFiles& w) const {
	if(w.is_binary()){
		w<<Vector<double>(1,mu_);
		w.add_to_header()->title(RST::math("\\mu")+"="+my::tostring(mu_),'<');
		w.add_to_header()->add(system_info_.get());
	} else { w<<mu_<<" "; }
}

Matrix<double> Squarek2MuPi::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 1.0;
	tmp(1,0) =-1.0;
	tmp(0,1) = 1.0;
	tmp(1,1) = 1.0;
	return tmp;
}

unsigned int Squarek2MuPi::unit_cell_index(Vector<double> const& x) const {
	if(my::are_equal(x(0),0.0,eq_prec_,eq_prec_) && my::are_equal(x(1),0.0,eq_prec_,eq_prec_)){ return 0; }
	if(my::are_equal(x(0),0.5,eq_prec_,eq_prec_) && my::are_equal(x(1),0.5,eq_prec_,eq_prec_)){ return 1; }
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void Squarek2MuPi::display_results(){
	compute_H(0);
	std::string mu(my::tostring(mu_));
	draw_lattice(true,true,false,dir_nn_[2]*0.5+dir_nn_[1],"-d:mu "+mu,RST::math("\\mu")+"="+mu);
}

void Squarek2MuPi::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="square-k2mupi";
	display_results();

	//plot_band_structure();
}
/*}*/
