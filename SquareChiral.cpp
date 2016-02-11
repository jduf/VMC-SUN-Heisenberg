#include "SquareChiral.hpp"

SquareChiral::SquareChiral(System const& s, double const& phi):
	System(s),
	Square<std::complex<double> >(set_ab(ref_(3),N_/m_),N_/m_,"square-chiral"),
	phi_(phi*M_PI/N_)
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		system_info_.text("SquareChiral :");
		system_info_.text(" Each color has the same Hamiltonian.");
		system_info_.text(" There is a flux of "+RST::math(my::tostring(phi)+"\\pi/"+my::tostring(N_))+ "per square plaquette");

		filename_ += "-phi"+my::tostring(phi);
	}
}

/*{method needed for running*/
void SquareChiral::compute_H(){
	H_.set(n_,n_,0);

	double t(1.0);
	unsigned int s0(0);
	unsigned int s1(0);
	
	switch(spuc_){
		case 4:
			{
				for(unsigned int i(0);i<obs_[0].nlinks();i++){
					s0 = obs_[0](i,0);
					s1 = obs_[0](i,1);
					switch(obs_[0](i,5)){
						case 0: { H_(s0,s1) = (obs_[0](i,4)?bc_*t:t); }break;
						case 1: { H_(s0,s1) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)?2.0*phi_:0.0)); }break;
						case 2: { H_(s0,s1) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)?0:phi_)); }break;
						case 3: { H_(s0,s1) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)?0.0:phi_)); }break;
					}
				}
			}break;
		case 5:
			{
				for(unsigned int i(0);i<obs_[0].nlinks();i++){
					s0 = obs_[0](i,0);
					s1 = obs_[0](i,1);
					switch(obs_[0](i,5)){
						case 0: { H_(s0,s1) = (obs_[0](i,4)?bc_*t:t); }break;
						case 1: { H_(s0,s1) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)?-phi_:phi_)); }break;
						case 2: { H_(s0,s1) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)?phi_:0.0)); }break;
						case 3: { H_(s0,s1) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)?0.0:-phi_)); }break;
						case 4: { H_(s0,s1) = (obs_[0](i,4)?bc_*t:t); }break;
					}
				}
			}break;
	}
	H_ += H_.conjugate_transpose();
}

void SquareChiral::create(){
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

void SquareChiral::save_param(IOFiles& w) const {
	std::string s("phi=("+my::tostring(phi_/(M_PI*m_))+")");
	Vector<double> param(1,phi_);

	w.add_header()->title(s,'<');
	w<<param;
	GenericSystem<std::complex<double> >::save_param(w);
}

Matrix<double> SquareChiral::set_ab(unsigned int const& ref3, unsigned int const& k) const {
	Matrix<double> tmp;
	switch(k){
		case 4:
			{
				tmp.set(2,2);
				tmp(0,0) = 2.0;
				tmp(1,0) = 0.0;
				tmp(0,1) = 0.0;
				tmp(1,1) = 2.0;
			}break;
		case 5:
			{
				tmp.set(2,2);
				if(ref3==2){
					tmp(0,0) = 2.0;
					tmp(1,0) = 1.0;
					tmp(0,1) =-1.0;
					tmp(1,1) = 2.0;
				} else {
					tmp(0,0) = 2.0;
					tmp(1,0) =-1.0;
					tmp(0,1) = 1.0;
					tmp(1,1) = 2.0;
				}
			}break;
	}
	return tmp;
}

unsigned int SquareChiral::match_pos_in_ab(Vector<double> const& x) const {
	Vector<double> match(2,0);
	switch(N_/m_){
		case 4:
			{
				if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 0; }
				match(0) = 0.5;
				if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 1; }
				match(0) = 0.0;
				match(1) = 0.5;
				if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 2; }
				match(0) = 0.5;
				match(1) = 0.5;
				if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 3; }
			}
		case 5:
			{
				if(ref_(3)==2){
					if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 0; }
					match(0) = 0.2;
					match(1) = 0.4;
					if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 1; }
					match(0) = 0.6;
					match(1) = 0.2;
					if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 2; }
					match(0) = 0.4;
					match(1) = 0.8;
					if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 3; }
					match(0) = 0.8;
					match(1) = 0.6;
					if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 4; }
				} else {
					if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 0; }
					match(0) = 0.4;
					match(1) = 0.2;
					if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 1; }
					match(0) = 0.8;
					match(1) = 0.4;
					if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 2; }
					match(0) = 0.2;
					match(1) = 0.6;
					if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 3; }
					match(0) = 0.6;
					match(1) = 0.8;
					if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 4; }
				}
			}
	}
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return 5;
}
/*}*/

/*{method needed for checking*/
void SquareChiral::display_results(){
	compute_H();

	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	std::string arrow("-");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(info_+path_+dir_,filename_);
	ps.begin(-20,-20,20,20,filename_);
	ps.polygon(cluster_vertex_,"linecolor=green");
	ps.polygon(draw_unit_cell(),"linecolor=black");

	std::complex<double> t;
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
				ps.put(xy1(0)+0.2,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
			} else { linestyle = "solid"; }

			if(t.real()>0){ color = "blue"; }
			else          { color = "red"; }

			if(my::are_equal(t.imag(),0.0)){
				arrow = "-";
			} else {
				if(t.imag()>0){ arrow = "->"; }
				else          { arrow = "<-"; }
			}

			ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}
		if(i%2){ ps.put(xy0(0)+0.10,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}"); }
	}
	ps.line("-",boundary_vertex_[0](0),boundary_vertex_[0](1),boundary_vertex_[1](0),boundary_vertex_[1](1),"linecolor=yellow");
	ps.line("-",boundary_vertex_[3](0),boundary_vertex_[3](1),boundary_vertex_[0](0),boundary_vertex_[0](1),"linecolor=yellow");
	ps.end(true,true,true);
}

void SquareChiral::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="square-chiral";
	display_results();

	//plot_band_structure();
	std::cout<<system_info_<<std::endl;
}
/*}*/
