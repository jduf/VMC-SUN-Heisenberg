#include "SquareFreeHopping.hpp"

SquareFreeHopping::SquareFreeHopping(System const& s, Vector<double> const& t):
	System(s),
	Square<double>(set_ab(),2,"square-free-real"),
	t_(t)
{
	if(status_==2){
		init_fermionic();

		system_info_.text("FreeReal : all colors experience the same Hamiltonian");
	}
}

/*{method needed for running*/
void SquareFreeHopping::compute_H(){
	H_.set(n_,n_,0);
	unsigned int k(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_:1)*t_(k);
		k = (k+1)%t_.size();
	}
	H_ += H_.transpose();
}

void SquareFreeHopping::create(){
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

Matrix<double> SquareFreeHopping::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 2;
	tmp(1,0) = 0;
	tmp(0,1) = 0;
	tmp(1,1) = 1;
	return tmp;
}

unsigned int SquareFreeHopping::match_pos_in_ab(Vector<double> const& x) const {
	Vector<double> match(2,0);
	if(my::are_equal(x,match)){ return 0; }
	match(0) = 0.5;
	match(1) = 0;
	if(my::are_equal(x,match)){ return 1; }
	return 2;
}
/*}*/

/*{method needed for checking*/
void SquareFreeHopping::display_results(){
	compute_H();
	std::string color("black");
	std::string linestyle("solid");
	std::string arrow("-");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(info_+path_+dir_,filename_);
	ps.begin(-9,-10,16,10,filename_);
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
	for(unsigned int i(0);i<polygon.row();i++){
		polygon(i,0) -= 0.2;
		polygon(i,1) -= 0.1;
	}
	ps.polygon(polygon,"linecolor=black");

	std::complex<double> t;
	unsigned int s0;
	unsigned int s1;
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		xy0 = x_[s0];

		s1 = obs_[0](i,1);
		xy1 = x_[s1];

		if((xy0-xy1).norm_squared()<1.1){ linestyle = "solid"; }
		else {
			linestyle = "dashed";
			xy1 = (xy0+dir_nn_[obs_[0](i,3)]).chop();
			ps.put(xy1(0)-0.20,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
		}

		t = H_(s0,s1);
		if(i%2){
			ps.put(xy0(0)-0.20,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}");
			if(my::real(H_(s0,s0))){ ps.circle(xy0,t.real(),"linecolor=magenta,fillstyle=solid,fillcolor=magenta"); }
		}

		if(std::abs(t)>1e-4){
			if(t.real()<0){ color = "red"; }
			else { color = "blue"; }

			arrow = "-";
			if(std::arg(t)>0){ arrow = "-"+std::string(std::arg(t)/(2*M_PI*m_/N_),'>'); }
			if(std::arg(t)<0){ arrow = std::string(-std::arg(t)/(2*M_PI*m_/N_),'<')+"-"; }
			ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+my::tostring(std::abs(t))+"pt,linecolor="+color+",linestyle="+linestyle);
		}
	}
	ps.end(true,true,true);
}

void SquareFreeHopping::check(){
	//unsigned int c(0);
	//unsigned int a(M_(c)-1);
	//unsigned int b(M_(c)-1);
	//Vector<double> eval;
	//do{ b++; } while (b+1<n_ && my::are_equal(eval(b),eval(b-1)));
	//if(b!=M_(c)){ while(a>0 && my::are_equal(eval(a-1),eval(a))){ a--; } }
	//std::cout<<a<<" "<<b<<std::endl;
	std::cout<<t_<<std::endl;
	compute_H();
	for(unsigned int i(0);i<n_;i++){
		for(unsigned int j(i);j<n_;j++){
			if(H_(i,j)!=0){ std::cout<<i<<" "<<j<<" "<<H_(i,j)<<std::endl; }
		}
	}
	//plot_band_structure();
	status_++;
	std::cout<<H_<<std::endl;
}
/*}*/
