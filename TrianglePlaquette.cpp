#include "TrianglePlaquette.hpp"

TrianglePlaquette::TrianglePlaquette(System const& s, double const& t):
	System(s),
	Triangle<double>(set_ab(),3,"triangle-plaquette"),
	t_(t)
{
	if(status_==2){
		init_fermionic();

		system_info_.text("Plaquette : all colors experience the same Hamiltonian");
		filename_ += "-t"+std::string(t_>=0?"+":"")+my::tostring(t_);
	}
}

/*{method needed for running*/
void TrianglePlaquette::compute_H(){
	H_.set(n_,n_,0);

	double t(-1.0);
	unsigned int s0(0);
	unsigned int s1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		switch(get_site_in_ab(s0)){
			case 0:
				{
					if(obs_[0](i,3)==2){ H_(s0,s1) = (obs_[0](i,4)?bc_*t_:t_); } 
					else { H_(s0,s1) = (obs_[0](i,4)?bc_*t:t); }
				}break;
			case 1:
				{
					if(obs_[0](i,3)!=2){ H_(s0,s1) = (obs_[0](i,4)?bc_*t_:t_); } 
					else { H_(s0,s1) = (obs_[0](i,4)?bc_*t:t); }
				}break;
			case 2:
				{ H_(s0,s1) = (obs_[0](i,4)?bc_*t_:t_); }break;
		}
	}
	H_ += H_.transpose();
}

void TrianglePlaquette::create(){
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

void TrianglePlaquette::save_param(IOFiles& w) const {
	std::string s("t=("+my::tostring(t_)+")");
	Vector<double> param(1,t_);

	w.add_header()->title(s,'<');
	w<<param;
	GenericSystem<double>::save_param(w);
}

unsigned int TrianglePlaquette::match_pos_in_ab(Vector<double> const& x) const {
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

Matrix<double> TrianglePlaquette::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 1.5;
	tmp(1,0) =-sqrt(3.0)/2;
	tmp(0,1) = 1.5;
	tmp(1,1) = sqrt(3.0)/2;
	return tmp;
}
/*}*/

/*{method needed for checking*/
void TrianglePlaquette::display_results(){
	compute_H();

	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1mm");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(info_+path_+dir_,filename_);
	ps.begin(-20,-20,20,20,filename_);
	ps.polygon(lattice_corners_,"linecolor=green");

	double x_shift(-(ab_(0,0)+ab_(0,1))/2);
	double y_shift(0.0);
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
				xy1 = xy0;
				xy1+= dir_nn_[obs_[0](i,3)];
				xy1 = xy1.chop();
				ps.put(xy1(0)+0.2,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
			} else { linestyle = "solid"; }

			if(t>0){ color = "blue"; }
			else   { color = "red"; }
			linewidth = my::tostring(std::abs(t))+"mm";
			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}
		if(i%3==2){ ps.put(xy0(0)+0.2,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}"); }
	}
	ps.line("-",boundary_[0](0),boundary_[0](1),boundary_[1](0),boundary_[1](1),"linecolor=yellow");
	ps.line("-",boundary_[3](0),boundary_[3](1),boundary_[0](0),boundary_[0](1),"linecolor=yellow");
	ps.end(true,true,true);
}

void TrianglePlaquette::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="triangle-plaquette";
	display_results();
	plot_band_structure();
}
/*}*/
