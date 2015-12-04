#include "SquareFreeFlux.hpp"

SquareFreeFlux::SquareFreeFlux(System const& s, Vector<double> const& phi):
	System(s),
	Square<std::complex<double> >((N_%m_?0:N_/m_),0,1,"square-free-flux"),
	phi_(phi)
{
	if(status_==2){
		init_fermionic();

		system_info_.text("FreeComplex : all colors experience the same Hamiltonian");
	}
}

/*{method needed for running*/
void SquareFreeFlux::compute_H(){
	H_.set(n_,n_,0);
	unsigned int s(0);
	unsigned int s0(0);
	unsigned int s1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		s = get_site_in_ab(s0);
		H_(s0,s1) = std::polar(obs_[0](i,3)*1.0,phi_(s));
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
	std::string linewidth("1pt");
	std::string arrow("-");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	std::complex<double> t;
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
				if(obs_[0](i,3)){ xy1(1) += 1.0; }
				else { xy1(0) += 1.0; }
				xy1 = xy1.chop();
				ps.put(xy1(0)-0.20,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
			} else { linestyle = "solid"; }

			if(t.real()>0){ color = "blue"; }
			else { color = "red"; }

			arrow = "-";
			if(std::arg(t)>0){ arrow = "-"+std::string(std::arg(t)/(2*M_PI*m_/N_),'>'); }
			if(std::arg(t)<0){ arrow = std::string(-std::arg(t)/(2*M_PI*m_/N_),'<')+"-"; }

			linewidth = my::tostring(std::abs(t))+"mm";
			ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}
		if(i%2){
			ps.put(xy0(0)-0.20,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}");
			if(my::real(H_(s0,s0))){ ps.circle(xy0,t.real(),"linecolor=magenta,fillstyle=solid,fillcolor=magenta"); }
		}
	}

	if(obs_.size()>1){
		double y_shift(4);
		double lr_corr;
		double rescale(obs_[1].nlinks()?0.75/obs_[1][0].get_x():0);
		for(unsigned int i(0);i<obs_[1].nlinks();i++){
			lr_corr = obs_[1][i].get_x()*rescale;
			if(std::abs(lr_corr)>1e-4){
				xy1 = get_pos_in_lattice(i);
				set_pos_LxLy(xy1);
				xy1 = (LxLy_*xy1).chop();
				xy1(1) -= 2*y_shift;

				if(i){
					if(lr_corr<0){ color = "red"; }
					else { color = "blue"; }
				} else { color = "black"; }

				ps.circle(xy1,std::abs(lr_corr),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
			}
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
