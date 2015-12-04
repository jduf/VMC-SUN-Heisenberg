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
	double t(1.0);

	unsigned int s0(0);
	unsigned int s1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		H_(s0,s1) = obs_[0](i,4)*t;
	}
	H_ += H_.transpose();
}

void TriangleFermi::create(){
	compute_H();
	diagonalize(false);
	for(unsigned int c(0);c<N_;c++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_(c);j++){
				EVec_[c](i,j) = H_(i,j);
			}
		}
	}
}

Matrix<double> TriangleFermi::set_ab(){
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
	std::string linewidth("1pt");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	double t;
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

	unsigned int s0;
	unsigned int s1;
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		xy0 = get_pos_in_lattice(s0);

		s1 = obs_[0](i,1);
		xy1 = get_pos_in_lattice(s1);

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
		ps.put(xy0(0)-0.20,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}");
	}
	ps.end(true,true,true);
}

void TriangleFermi::check(){
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
	filename_ ="triangle-fermi";
	display_results();
}
/*}*/
