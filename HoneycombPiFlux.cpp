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
	Matrix<int> nb;
	unsigned int s(0);
	for(unsigned int i(0);i<n_;i+=2){
		s = get_site_in_ab(i);
		nb = get_neighbourg(i);
		switch(s){
			case 0:
				{
					H_(i,nb(0,0))= nb(0,1)*th;
					H_(i,nb(1,0))= nb(1,1)*th;
					H_(i,nb(2,0))= nb(2,1)*th;
				}break;
			case 2:
				{
					H_(i,nb(0,0))= nb(0,1)*th;
					H_(i,nb(1,0))= nb(1,1)*th;
					H_(i,nb(2,0))= nb(2,1)*td;
				}break;
		}
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

unsigned int HoneycombPiFlux::match_pos_in_ab(Vector<double> const& x) const{
	Vector<double> match(2,0);
	if(my::are_equal(x,match)){ return 0; }
	match(0) = 1.0/3.0;
	if(my::are_equal(x,match)){ return 1; }
	match(0) = 0.5;
	match(1) = 0.5;
	if(my::are_equal(x,match)){ return 2; }
	match(0) = 1.0/3.0;
	if(my::are_equal(x,match)){ return 3; }
	return 4;
}

Matrix<double> HoneycombPiFlux::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 3.0;
	tmp(1,0) = 0;
	tmp(0,1) = 0.0;
	tmp(1,1) = sqrt(3.0);
	return tmp;
}
/*}*/

/*{method needed for checking*/
void HoneycombPiFlux::display_results(){
	compute_H();

	Matrix<int> nb;
	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	double t;
	PSTricks ps(info_+path_+dir_,filename_);
	ps.begin(-20,-10,20,10,filename_);

	Matrix<double> polygon(4,2);
	polygon(0,0)=0;
	polygon(0,1)=0;
	polygon(1,0)=LxLy_(0,0);
	polygon(1,1)=LxLy_(1,0);
	polygon(2,0)=LxLy_(0,0)+LxLy_(0,1);
	polygon(2,1)=LxLy_(1,0)+LxLy_(1,1);
	polygon(3,0)=LxLy_(0,1);
	polygon(3,1)=LxLy_(1,1);
	ps.polygon(polygon,"linecolor=green");

	polygon(0,0)=0;
	polygon(0,1)=0;
	polygon(1,0)=ab_(0,0);
	polygon(1,1)=ab_(1,0);
	polygon(2,0)=ab_(0,0)+ab_(0,1);
	polygon(2,1)=ab_(1,0)+ab_(1,1);
	polygon(3,0)=ab_(0,1);
	polygon(3,1)=ab_(1,1);
	ps.polygon(polygon,"linecolor=black");

	for(unsigned int i(0);i<n_;i++){
		xy0 = get_pos_in_lattice(i);
		set_pos_LxLy(xy0);
		xy0 = (LxLy_*xy0).chop();
		ps.put(xy0(0)-0.20,xy0(1)+0.15,my::tostring(i));

		if(!(i%2)){
			nb = get_neighbourg(i);

			t = H_(i,nb(0,0));
			if(std::abs(t)>1e-5){
				xy1 = get_pos_in_lattice(nb(0,0));
				set_pos_LxLy(xy1);
				xy1 = LxLy_*xy1;
				if((xy0-xy1).norm_squared()>1.0001){
					linestyle = "dashed"; 
					xy1 = xy0;
					xy1(0) += dir_nn_(0,0);
					xy1(1) += dir_nn_(0,1);
					ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(nb(0,0)));
				} else { linestyle = "solid";  }
				xy1 = xy1.chop();

				if(t>0){ color = "blue";}
				else { color = "red"; }
				linewidth = my::tostring(std::abs(t))+"mm";
				/*(+x)-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
			}

			t = H_(i,nb(1,0));
			if(std::abs(t)>1e-5){
				xy1 = get_pos_in_lattice(nb(1,0));
				set_pos_LxLy(xy1);
				xy1 = LxLy_*xy1;
				if((xy0-xy1).norm_squared()>1.0001){
					linestyle = "dashed"; 
					xy1 = xy0;
					xy1(0) += dir_nn_(1,0);
					xy1(1) += dir_nn_(1,1);
					ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(nb(1,0)));
				} else { linestyle = "solid";  }
				xy1 = xy1.chop();

				if(t>0){ color = "blue";}
				else { color = "red"; }
				linewidth = my::tostring(std::abs(t))+"mm";
				/*(+y)-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
			}

			t = H_(i,nb(2,0));
			if(std::abs(t)>1e-5){
				xy1 = get_pos_in_lattice(nb(2,0));
				set_pos_LxLy(xy1);
				xy1 = LxLy_*xy1;
				if((xy0-xy1).norm_squared()>1.0001){
					linestyle = "dashed"; 
					xy1 = xy0;
					xy1(0) += dir_nn_(2,0);
					xy1(1) += dir_nn_(2,1);
					ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(nb(2,0)));
				} else { linestyle = "solid";  }
				xy1 = xy1.chop();

				if(t>0){ color = "blue";}
				else { color = "red"; }
				linewidth = my::tostring(std::abs(t))+"mm";
				/*(-y)-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
			}
		}
	}
	ps.end(true,true,true);
}

void HoneycombPiFlux::check(){
	info_ ="";
	path_ ="";
	dir_ ="./";
	filename_ ="honeycomb-piflux";
	display_results();
}
/*}*/
