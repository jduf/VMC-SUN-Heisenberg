#include "HoneycombPiFlux.hpp"

HoneycombPiFlux::HoneycombPiFlux(System const& s):
	System(s),
	Honeycomb<double>(set_ab(),4,"honeycomb-piflux")
{
	if(status_==2){
		init_fermionic();

		system_info_.item("pi-flux configuration");
		system_info_.item("4 sites per unit cell");
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
		if(get_site_in_ab(s0)==2 && get_site_in_ab(s1)==1 && obs_[0](i,3)==2){ H_(s0,s1) = (obs_[0](i,4)?bc_:1)*td; }
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

unsigned int HoneycombPiFlux::match_pos_in_ab(Vector<double> const& x) const {
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

	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	double t;
	PSTricks ps(info_+path_+dir_,filename_);
	ps.begin(-20,-20,20,20,filename_);
	ps.polygon(lattice_corners_,"linecolor=green");

	Matrix<double> polygon(4,2);
	polygon(0,0)=0;
	polygon(0,1)=0;
	polygon(1,0)=ab_(0,0);
	polygon(1,1)=ab_(1,0);
	polygon(2,0)=ab_(0,0)+ab_(0,1);
	polygon(2,1)=ab_(1,0)+ab_(1,1);
	polygon(3,0)=ab_(0,1);
	polygon(3,1)=ab_(1,1);
	ps.polygon(polygon,"linecolor=black");

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
				ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(s1));
			} else { 
				linestyle = "solid";  
				if(s0<s1){
					ps.put(xy0(0)-0.20,xy0(1)+0.15,my::tostring(s0)); 
					ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(s1)); 
				}
			}

			if(t>0){ color = "blue";}
			else   { color = "red"; }
			linewidth = my::tostring(std::abs(t))+"mm";
			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}
	}
	ps.end(true,true,true);
}

void HoneycombPiFlux::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="honeycomb-piflux";
	display_results();
}
/*}*/
