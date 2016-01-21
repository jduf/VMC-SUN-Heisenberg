#include "TriangleFree.hpp"

TriangleFree::TriangleFree(System const& s, Vector<double> const& t, Vector<double> const& mu):
	System(s),
	Triangle<double>(set_ab(),3,"triangle-free"),
	t_(t),
	mu_(mu)
{
	if(status_==2){
		init_fermionic();

		system_info_.text("Free : all colors experience the same Hamiltonian");
	}
}

/*{method needed for running*/
void TriangleFree::compute_H(){
	H_.set(n_,n_,0);

	unsigned int s0(0);
	unsigned int s1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);

		switch(get_site_in_ab(s0)){
			case 0:
				{
					if(obs_[0](i,3)==2){ H_(s0,s1) = (obs_[0](i,4)?bc_:1)*t_(0); } 
					else { H_(s0,s1) = (obs_[0](i,4)?bc_:1)*t_(1); }
					H_(s0,s0) = mu_(0);
				}break;
			case 1:
				{
					if(obs_[0](i,3)!=2){ H_(s0,s1) = (obs_[0](i,4)?bc_:1)*t_(0); } 
					else { H_(s0,s1) = (obs_[0](i,4)?bc_:1)*t_(1); }
					H_(s0,s0) = mu_(1);
				}break;
			case 2:
				{ 
					H_(s0,s1) = (obs_[0](i,4)?bc_:1)*t_(0); 
					H_(s0,s0) = mu_(2);
				}break;
		}
	}
	H_ += H_.transpose();
}

void TriangleFree::create(){
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

unsigned int TriangleFree::match_pos_in_ab(Vector<double> const& x) const {
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

Matrix<double> TriangleFree::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 1.5;
	tmp(1,0) =-sqrt(3.0)/2;
	tmp(0,1) = 1.5;
	tmp(1,1) = sqrt(3.0)/2;
	return tmp;
}
/*}*/

/*{method needed for checking*/
void TriangleFree::display_results(){
	compute_H();

	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(info_+path_+dir_,filename_);
	ps.begin(-2,-20,20,20,filename_);

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

	//double x_shift((LxLy_(0,0)+LxLy_(0,1)-ab_(0,0)-ab_(0,1))/2);
	//double y_shift((ab_(1,0)+ab_(1,1)/2)/2);
	//double x_shift(-ab_(0,0)/4);
	//double y_shift(0.0);
	//polygon(0,0)+=x_shift;
	//polygon(0,1)+=y_shift;
	//polygon(1,0)+=x_shift;
	//polygon(1,1)+=y_shift;
	//polygon(2,0)+=x_shift;
	//polygon(2,1)+=y_shift;
	//polygon(3,0)+=x_shift;
	//polygon(3,1)+=y_shift;
	ps.polygon(polygon,"linecolor=black");

	double t;
	double mu;
	unsigned int s0;
	unsigned int s1;
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		xy0 = get_pos_in_lattice(s0);

		s1 = obs_[0](i,1);
		xy1 = get_pos_in_lattice(s1);

		//if((my::in_polygon(polygon.row(),polygon.ptr(),polygon.ptr()+polygon.row(),xy0(0),xy0(1)) || my::in_polygon(polygon.row(),polygon.ptr(),polygon.ptr()+polygon.row(),xy1(0),xy1(1))) ){
		t = H_(s0,s1);
		if(std::abs(t)>1e-4){
			if((xy0-xy1).norm_squared()>1.0001){
				linestyle = "dashed";
				xy1 = xy0;
				xy1(0) += dir_nn_(obs_[0](i,3),0);
				xy1(1) += dir_nn_(obs_[0](i,3),1);
				xy1 = xy1.chop();
				ps.put(xy1(0)-0.20,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
			} else { linestyle = "solid"; }

			if(t>0){ color = "blue"; }
			else { color = "red"; }
			linewidth = my::tostring(std::abs(t))+"mm";
			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}

		mu = H_(s0,s0);
		if(std::abs(mu)>1e-4){
			if(mu<0){ color = "magenta"; }
			else { color = "cyan"; }
			ps.circle(xy0,std::abs(mu),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
		}
		//}
		ps.put(xy0(0)-0.20,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}");
	}
	ps.end(true,true,true);
}

void TriangleFree::check(){
	//Matrix<int> nb;
	//std::cout<<"######################"<<std::endl;
	//for(unsigned int i(0);i<n_;i++){
	//nb = get_neighbourg(i);
	//std::cout<<"i="<<i<<std::endl;
	//std::cout<<nb<<std::endl;
	//}
	//std::cout<<"######################"<<std::endl;
	//nb = get_neighbourg(3);
	//std::cout<<nb<<std::endl;
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="triangle-free";
	display_results();
}
/*}*/
