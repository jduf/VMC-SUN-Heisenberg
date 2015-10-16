#include "SquareACSL.hpp"

SquareACSL::SquareACSL(System const& s, Vector<double> const& t):
	System(s),
	Square<std::complex<double> >((N_%m_?0:N_/m_),N_/m_,0,"square-acsl"),
	t_(t)
{
	if(2*spuc_ != t_.size()){
		status_++;
		std::cerr<<__PRETTY_FUNCTION__<<" : t has a wrong size"<<std::endl; 
	}
	if(status_==2){
		init_fermionic();

		system_info_.text("ACSL : all colors experience the same Hamiltonian");
	}
}

/*{method needed for running*/
void SquareACSL::compute_H(){
	H_.set(n_,n_,0);
	Matrix<int> nb;
	unsigned int s(0);
	double phi(2*M_PI*m_/N_);
	for(unsigned int i(0);i<n_;i++){
		s = get_site_in_ab(i);
		nb = get_neighbourg(i);
		if(s){
			H_(i,nb(0,0)) = t_(2*s-1)*nb(0,1);
			H_(i,nb(1,0)) = std::polar(t_(2*s)*nb(1,1),s*phi);
		} else {
			H_(i,nb(0,0)) = nb(0,1);
			H_(i,nb(1,0)) = std::polar(t_(0)*nb(1,1),s*phi);
		}
	}
	H_ += H_.trans_conj();
}

void SquareACSL::create(){
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

unsigned int SquareACSL::match_pos_in_ab(Vector<double> const& x) const{
	Vector<double> match(2,0);
	for(unsigned int i(0);i<spuc_;i++){
		match(0) = 1.0/spuc_*i;
		if(my::are_equal(x,match)){ return i; }
	}
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void SquareACSL::lattice(std::string const& path, std::string const& filename){
	compute_H();
	std::string color("black");
	std::string linestyle("solid");
	std::string arrow("-");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(path,filename);
	ps.begin(-9,-10,16,10,filename_);
	std::complex<double> t;
	unsigned int s0;
	unsigned int s1;
	double y_shift(4);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		xy0 = get_pos_in_lattice(s0);
		set_pos_LxLy(xy0);
		xy0 = (LxLy_*xy0).chop();

		s1 = obs_[0](i,1);
		xy1 = get_pos_in_lattice(s1);
		set_pos_LxLy(xy1);
		xy1 = (LxLy_*xy1).chop();

		if((xy0-xy1).norm_squared()<1.1){ linestyle = "solid"; }
		else {
			linestyle = "dashed";
			if(i%2 && xy1(1)<xy0(1)){
				xy1(0) = xy0(0);
				xy1(1) = xy0(1)+1.0;
			}
			if(!(i%2) && xy1(0)<xy0(0)){
				xy1(0) = xy0(0)+1.0;
				xy1(1) = xy0(1);
			}
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
			
			xy0 = xy0.chop();
			xy1 = xy1.chop();
			ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+my::tostring(std::abs(t))+"pt,linecolor="+color+",linestyle="+linestyle);
		}
	}

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
	for(unsigned int i(0);i<polygon.row();i++){
		polygon(i,0) -= 0.2;
		polygon(i,1) -= 0.1;
	}
	ps.polygon(polygon,"linecolor=black");
	ps.end(true,true,true);
}

void SquareACSL::check(){
	lattice("./","lattice");
}
/*}*/
