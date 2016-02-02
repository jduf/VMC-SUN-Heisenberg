#include "SquareFreeFlux.hpp"

SquareFreeFlux::SquareFreeFlux(System const& s, Vector<double> const& phi):
	System(s),
	Square<std::complex<double> >(set_ab(),4,"square-freeflux"),
	phi_(phi)
{
	if(status_==2){
		init_lattice();
		init_fermionic();

		system_info_.text("FreeComplex : all colors experience the same Hamiltonian");
	}
}

/*{method needed for running*/
void SquareFreeFlux::compute_H(){
	H_.set(n_,n_,0);
	unsigned int s0(0);
	unsigned int s1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		if(obs_[0](i,3)==1){ H_(s0,s1) = my::chop(std::polar(double(obs_[0](i,4)?bc_:1),phi_(obs_[0](i,5)))); }
		else { H_(s0,s1) = (obs_[0](i,4)?bc_:1); }
	}
	H_ += H_.conjugate_transpose();
}

void SquareFreeFlux::create(){
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

Matrix<double> SquareFreeFlux::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 2;
	tmp(1,0) = 0;
	tmp(0,1) = 0;
	tmp(1,1) = 2;
	return tmp;
}

unsigned int SquareFreeFlux::match_pos_in_ab(Vector<double> const& x) const {
	Vector<double> match(2,0);
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 0; }
	match(0) = 0.5;
	match(1) = 0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 1; }
	match(0) = 0;
	match(1) = 0.5;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 2; }
	match(0) = 0.5;
	match(1) = 0.5;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 3; }
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return 4;
}
/*}*/

/*{method needed for checking*/
void SquareFreeFlux::display_results(){
	compute_H();

	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1mm");
	std::string arrow("-");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	std::complex<double> t;
	PSTricks ps(info_+path_+dir_,filename_);
	ps.begin(-2,-20,20,20,filename_);
	ps.polygon(lattice_corners_,"linecolor=green");

	Matrix<double> polygon(4,2);
	polygon(0,0) = 0;
	polygon(0,1) = 0;
	polygon(1,0) = ab_(0,0);
	polygon(1,1) = ab_(1,0);
	polygon(2,0) = ab_(0,0)+ab_(0,1);
	polygon(2,1) = ab_(1,0)+ab_(1,1);
	polygon(3,0) = ab_(0,1);
	polygon(3,1) = ab_(1,1);
	polygon -= 0.1;
	ps.polygon(polygon,"linecolor=black");

	unsigned int s0;
	unsigned int s1;
	double phi(0);
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
				ps.put(xy1(0)-0.20,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
			} else { linestyle = "solid"; }

			if(t.real()>0){ color = "blue"; }
			else { color = "red"; }

			arrow = "-";
			if(std::arg(t)>0){ arrow = "-"+std::string(std::arg(t)/(2*M_PI),'>'); }
			if(std::arg(t)<0){ arrow = std::string(-std::arg(t)/(2*M_PI),'<')+"-"; }

			linewidth = my::tostring(std::abs(t))+"mm";
			ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}
		if(i%2){
			ps.put(xy0(0)+0.10,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}");
			if(my::real(H_(s0,s0))){ ps.circle(xy0,t.real(),"linecolor=magenta,fillstyle=solid,fillcolor=magenta"); }
			ps.put(my::chop((2*xy0(0)-1)/2.0),my::chop((xy0(1)+xy1(1))/2.0),"\\tiny{"+my::tostring((std::arg(t)-phi)/(2*M_PI))+"}");
			phi =  std::arg(t);
		}
	}
	ps.end(true,true,true);
}

void SquareFreeFlux::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="square-freeflux";
	display_results();
}
/*}*/
