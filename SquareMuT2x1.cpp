#include "SquareMuT2x1.hpp"

SquareMuT2x1::SquareMuT2x1(System const& s, double const& mu, Vector<double> const& t):
	System(s),
	Square<double>(set_ab(),2,"square-mu-T2x1"),
	mu_(mu),
	t_(t)
{
	if(N_/m_==2 && t_.size() == 4){
		if(status_==3){ init_lattice(); }
		if(status_==2){
			init_fermionic();
			same_wf_ = false;

			system_info_.text("SquareMuT2x1 :");
			system_info_.item("Each color has a different Hamiltonian.");
			system_info_.item("Two sites per unit cell (2x1)");

			filename_ += "-mu"+std::string(mu_>=0?"+":"")+my::tostring(mu_);
			filename_ += "-t";
			for(unsigned int i(0);i<t_.size();i++){
				filename_ += ((t_(i)>=0)?"+":"")+my::tostring(t_(i));
			}
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : N/m != 2 or t_.size()!=4"<<std::endl; }
}

/*{method needed for running*/
void SquareMuT2x1::compute_H(unsigned int const& c){
	H_.set(n_,n_,0);

	unsigned int s0(0);
	unsigned int s1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		H_(s0,s1) = (obs_[0](i,4)?bc_*t_(obs_[0](i,2)):t_(obs_[0](i,2)));
		if(obs_[0](i,3) && (unsigned int)(obs_[0](i,5))==c%2){ H_(s0,s0) = mu_/2; }
	}
	H_ += H_.transpose();
}

void SquareMuT2x1::create(){
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

void SquareMuT2x1::save_param(IOFiles& w) const {
	if(w.is_binary()){
		Vector<double> param(t_.size()+1);
		param(0) = mu_;
		for(unsigned int i(0);i<t_.size();i++){ param(i+1) = t_(i); }
		w<<param;
		w.add_to_header()->title(RST::math("\\mu")+"="+my::tostring(mu_)+" t=("+my::tostring(t_)+")",'<');
		w.add_to_header()->add(system_info_.get());
	} else { w<<mu_<<" "<<t_; }
}

Matrix<double> SquareMuT2x1::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 1.0;
	tmp(1,0) = 1.0;
	tmp(0,1) = 1.0;
	tmp(1,1) =-1.0;
	return tmp;
}

unsigned int SquareMuT2x1::unit_cell_index(Vector<double> const& x) const {
	if(my::are_equal(x(0),0.0,eq_prec_,eq_prec_) && my::are_equal(x(1),0.0,eq_prec_,eq_prec_)){ return 0; }
	if(my::are_equal(x(0),0.5,eq_prec_,eq_prec_) && my::are_equal(x(1),0.5,eq_prec_,eq_prec_)){ return 1; }
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void SquareMuT2x1::display_results(){
	compute_H(0);
	std::string t(my::tostring(t_));
	std::string mu(my::tostring(mu_));
	draw_lattice(true,true,false,dir_nn_[2]*0.5+dir_nn_[3]*0.5,"-d:t "+t+" -d:mu "+mu,"t=("+t+") "+RST::math("\\mu")+"="+mu);
}

void SquareMuT2x1::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="square-mu-T2x1";
	display_results();

	//plot_band_structure();
}
/*}*/
