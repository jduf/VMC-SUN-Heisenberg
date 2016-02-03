#include "SquareMu.hpp"

SquareMu::SquareMu(System const& s, double const& mu):
	System(s),
	Square<double>(set_ab(),5,"square-mu"),
	mu_(mu)
{
	if(status_==2){
		init_lattice();
		init_fermionic();
		same_wf_ = false;

		system_info_.text("Mu : each color has a different Hamiltonian");
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
		compute_H(c);
		status_ = 2;
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

void SquareMu::save_param(IOFiles& w) const {
	std::string s("mu=("+my::tostring(mu_)+")");
	Vector<double> param(1,mu_);

	w.add_header()->title(s,'<');
	w<<param;
	GenericSystem<double>::save_param(w);
}

Matrix<double> SquareMu::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 2.0;
	tmp(1,0) =-1.0;
	tmp(0,1) = 1.0;
	tmp(1,1) = 2.0;
	return tmp;
}

unsigned int SquareMu::match_pos_in_ab(Vector<double> const& x) const {
	Vector<double> match(2,0);
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
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return 5;
}
/*}*/

/*{method needed for checking*/
void SquareMu::display_results(){
	compute_H(0);

	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1mm");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(info_+path_+dir_,filename_);
	ps.begin(-20,-20,20,20,filename_);
	ps.polygon(lattice_corners_,"linecolor=green");

	double x_shift(-(ab_(0,0)+ab_(0,1))/2.0);
	double y_shift(-(ab_(1,0)+ab_(1,1))/2.0);
	Matrix<double> polygon(4,2);
	polygon(0,0)=x_shift;
	polygon(0,1)=y_shift;
	polygon(1,0)=x_shift+ab_(0,0);
	polygon(1,1)=y_shift+ab_(1,0);
	polygon(2,0)=x_shift+ab_(0,0)+ab_(0,1);
	polygon(2,1)=y_shift+ab_(1,0)+ab_(1,1);
	polygon(3,0)=x_shift+ab_(0,1);
	polygon(3,1)=y_shift+ab_(1,1);
	ps.polygon(polygon,"linecolor=black");

	double t;
	double mu;
	unsigned int s0;
	unsigned int s1;
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		xy0 = x_[s0];

		s1 = obs_[0](i,1);
		xy1 = x_[s1];

		t = H_(s0,s1);
		if(std::abs(t)>1e-4){
			if((xy0-xy1).norm_squared()>1.0001){
				linestyle = "dashed";
				xy1 = (xy0+dir_nn_[obs_[0](i,3)]).chop();
				ps.put(xy1(0)+0.2,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
			} else { linestyle = "solid"; }

			if(t>0){ color = "blue"; }
			else   { color = "red"; }
			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}

		mu = H_(s0,s0);
		if(std::abs(mu)>1e-4){
			if(mu<0){ color = "green"; }
			else    { color = "cyan"; }
			ps.circle(xy0,std::abs(mu),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
		}

		if(i%2){ ps.put(xy0(0)+0.2,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}"); }
	}
	ps.line("-",boundary_[0](0),boundary_[0](1),boundary_[1](0),boundary_[1](1),"linecolor=yellow");
	ps.line("-",boundary_[3](0),boundary_[3](1),boundary_[0](0),boundary_[0](1),"linecolor=yellow");
	ps.end(true,true,true);
}

void SquareMu::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="square-mu";
	display_results();
	plot_band_structure();
}
/*}*/
