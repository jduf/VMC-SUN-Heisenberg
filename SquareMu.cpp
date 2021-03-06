#include "SquareMu.hpp"

SquareMu::SquareMu(System const& s, double const& mu):
	System(s),
	Square<double>(set_ab(ref_(3),N_/m_),N_/m_,"square-mu"),
	mu_(mu)
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();
		same_wf_ = false;

		system_info_.text("SquareMu :");
		system_info_.item("Each color has a different Hamiltonian.");

		filename_ += "-mu"+std::string(mu_>=0?"+":"")+my::tostring(mu_);
	}
}

/*{method needed for running*/
void SquareMu::compute_H(unsigned int const& c){
	H_.set(n_,n_,0);

	double t(-1.0);
	unsigned int s0(0);
	unsigned int s1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		H_(s0,s1) = (obs_[0](i,4)?bc_*t:t);
		if((unsigned int)(obs_[0](i,5))==c%spuc_){ H_(s0,s0) = mu_/2; }
	}
	H_ += H_.transpose();
}

void SquareMu::create(){
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

void SquareMu::save_param(IOFiles& w) const {
	if(w.is_binary()){
		w<<Vector<double>(1,mu_);
		w.add_to_header()->title(RST::math("\\mu")+"="+my::tostring(mu_),'<');
		w.add_to_header()->add(system_info_.get());
	} else { w<<mu_<<" "; }
}

Matrix<double> SquareMu::set_ab(unsigned int const& ref3, unsigned int const& k) const {
	Matrix<double> tmp(2,2);
	if(k==2){
		tmp(0,0) = 1.0;
		tmp(1,0) = 1.0;
		tmp(0,1) = 1.0;
		tmp(1,1) =-1.0;
	}
	if(k==3){
		tmp(0,0) = 3.0;
		tmp(1,0) = 0.0;
		tmp(0,1) = 1.0;
		tmp(1,1) = 1.0;
	}
	if(k==5){
		if(ref3==2){
			tmp(0,0) = 2.0;
			tmp(1,0) = 1.0;
			tmp(0,1) =-1.0;
			tmp(1,1) = 2.0;
		} else {
			tmp(0,0) = 2.0;
			tmp(1,0) =-1.0;
			tmp(0,1) = 1.0;
			tmp(1,1) = 2.0;
		}
	}
	return tmp;
}

unsigned int SquareMu::unit_cell_index(Vector<double> const& x) const {
	if(spuc_ == 2){
		if(my::are_equal(x(0),0.0,eq_prec_,eq_prec_) && my::are_equal(x(1),0.0,eq_prec_,eq_prec_)){ return 0; }
		if(my::are_equal(x(0),0.5,eq_prec_,eq_prec_) && my::are_equal(x(1),0.5,eq_prec_,eq_prec_)){ return 1; }
	}
	if(spuc_ == 3){
		if(my::are_equal(x(1),0.0,eq_prec_,eq_prec_)){
			if(my::are_equal(x(0),0.0,    eq_prec_,eq_prec_)){ return 0; }
			if(my::are_equal(x(0),1.0/3.0,eq_prec_,eq_prec_)){ return 1; }
			if(my::are_equal(x(0),2.0/3.0,eq_prec_,eq_prec_)){ return 2; }
		}
	}
	if(spuc_ == 5){
		Vector<double> match(2,0);
		if(ref_(3)==2){
			if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 0; }
			match(0) = 0.2;
			match(1) = 0.4;
			if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 1; }
			match(0) = 0.4;
			match(1) = 0.8;
			if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 2; }
			match(0) = 0.6;
			match(1) = 0.2;
			if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 3; }
			match(0) = 0.8;
			match(1) = 0.6;
			if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 4; }
		} else {
			if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 0; }
			match(0) = 0.4;
			match(1) = 0.2;
			if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 1; }
			match(0) = 0.8;
			match(1) = 0.4;
			if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 2; }
			match(0) = 0.2;
			match(1) = 0.6;
			if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 3; }
			match(0) = 0.6;
			match(1) = 0.8;
			if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 4; }
		}
	}
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void SquareMu::display_results(){
	compute_H(0);
	std::string mu(my::tostring(mu_));
	draw_lattice(true,true,false,(spuc_==3?dir_nn_[2]*0.25+dir_nn_[3]*0.5:dir_nn_[2]*0.5+dir_nn_[1]),"-d:mu "+mu,RST::math("\\mu")+"="+mu);
}

void SquareMu::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="square-mu";
	display_results();

	//plot_band_structure();
}
/*}*/
