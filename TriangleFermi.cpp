#include "TriangleFermi.hpp"

TriangleFermi::TriangleFermi(System const& s):
	System(s),
	Triangle<double>(set_ab(),1,"triangle-fermi")
{
	if(status_==2){
		init_fermionic();

		system_info_.text("Fermi : all colors experience the same Hamiltonian");
	}
}

/*{method needed for running*/
void TriangleFermi::compute_H(){
	H_.set(n_,n_,0);
	double t(-1.0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*t:t);
	}
	H_ += H_.transpose();
}

void TriangleFermi::create(){
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

Matrix<double> TriangleFermi::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 1;
	tmp(1,0) = 0;
	tmp(0,1) = 0.5;
	tmp(1,1) = sqrt(3.0)/2.0;
	return tmp;
}
/*}*/

/*{method needed for checking*/
void TriangleFermi::display_results(){
	compute_H();

	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1mm");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(info_+path_+dir_,filename_);
	ps.begin(-20,-20,20,20,filename_);
	ps.polygon(lattice_corners_,"linecolor=green");

	double x_shift(-(ab_(0,0)+ab_(0,1))/2.0);
	double y_shift(-ab_(1,1)/2.0);
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
				xy1 = (xy0+dir_nn_[obs_[0](i,3)]).chop();
				ps.put(xy1(0)-0.2,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
			} else { linestyle = "solid"; }

			if(t>0){ color = "blue"; }
			else   { color = "red"; }
			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}
		if(i%3==2){ ps.put(xy0(0)+0.2,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}"); }
	}
	ps.line("-",boundary_[0](0),boundary_[0](1),boundary_[1](0),boundary_[1](1),"linecolor=yellow");
	ps.line("-",boundary_[3](0),boundary_[3](1),boundary_[0](0),boundary_[0](1),"linecolor=yellow");
	ps.end(true,true,true);
}

void TriangleFermi::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="triangle-fermi";
	display_results();
	compute_H();
	plot_band_structure();
}
/*}*/
