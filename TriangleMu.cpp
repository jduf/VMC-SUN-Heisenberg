#include "TriangleMu.hpp"

TriangleMu::TriangleMu(System const& s, double const& mu):
	System(s),
	Triangle<double>(set_ab(),3,"triangle-mu"),
	mu_(mu)
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();
		same_wf_ = false;

		system_info_.text("TriangleMu :");
		system_info_.item("Each color has a different Hamiltonian.");

		filename_ += "-mu"+std::string(mu_>=0?"+":"")+my::tostring(mu_);
	}
}

/*{method needed for running*/
void TriangleMu::compute_H(unsigned int const& c){
	H_.set(n_,n_,0);

	double t(-1.0);
	unsigned int s0(0);
	unsigned int s1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		H_(s0,s1) = (obs_[0](i,4)?bc_*t:t);
		if((unsigned int)(obs_[0](i,5))==c%spuc_){ H_(s0,s0) = mu_/2.0; }
	}
	H_ += H_.transpose();
}

void TriangleMu::create(){
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
		}
	}
}

void TriangleMu::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("mu=("+my::tostring(mu_)+")");
		Vector<double> param(1,mu_);

		w.add_to_header()->title(s,'<');
		w<<param;
		w.add_to_header()->add(system_info_.get());
	} else { w<<mu_<<" "; }
}

Matrix<double> TriangleMu::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 1.5;
	tmp(1,0) =-sqrt(3.0)/2;
	tmp(0,1) = 1.5;
	tmp(1,1) = sqrt(3.0)/2;
	return tmp;
}

unsigned int TriangleMu::unit_cell_index(Vector<double> const& x) const {
	Vector<double> match(2,0);
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 0; }
	match(0) = 1.0/3.0;
	match(1) = 1.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 1; }
	match(0) = 2.0/3.0;
	match(1) = 2.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 2; }
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return 3;
}
/*}*/

/*{method needed for checking*/
void TriangleMu::display_results(){
	compute_H(0);
	std::string mu(my::tostring(mu_));
	draw_lattice(true,true,false,dir_nn_[3]*0.5,"-d:mu "+mu,RST::math("\\mu")+"="+mu);
}

void TriangleMu::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="triangle-mu";
	//display_results();

	compute_H(0);
	plot_band_structure();
}
/*}*/
