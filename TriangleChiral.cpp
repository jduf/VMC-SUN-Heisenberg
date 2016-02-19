#include "TriangleChiral.hpp"

TriangleChiral::TriangleChiral(System const& s, double const& phi):
	System(s),
	Triangle<std::complex<double> >(set_ab(),3,"triangle-chiral"),
	phi_(phi)
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		system_info_.text("TriangleChiral :");
		system_info_.item("Each color has the same Hamiltonian.");
		system_info_.item("There is a flux of "+RST::math(my::tostring(phi)+"\\pi/3") + "per triangular plaquette.");

		filename_ += "-phi"+my::tostring(phi_);
	}
}

/*{method needed for running*/
void TriangleChiral::compute_H(){
	H_.set(n_,n_,0);

	double t(-1.0);
	unsigned int s0(0);
	unsigned int s1(0);
	double phi(phi_*M_PI/3.0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);

		switch(obs_[0](i,3)){
			case 0:{ H_(s0,s1) = (obs_[0](i,4)?bc_*t:t); }break;
			case 1:{ H_(s0,s1) = std::polar((obs_[0](i,4)?bc_*t:t),phi*(2.0*obs_[0](i,5)+1)); }break;
			case 2:{ H_(s0,s1) = std::polar((obs_[0](i,4)?bc_*t:t),phi* 2.0*obs_[0](i,5)); }break;
		}
	}
	H_ += H_.conjugate_transpose();
}

void TriangleChiral::create(){
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

void TriangleChiral::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("phi=("+my::tostring(phi_)+")");
		Vector<double> param(1,phi_);

		w.add_header()->title(s,'<');
		w<<param;
		GenericSystem<std::complex<double> >::save_param(w);
	} else { w<<phi_<<" "; }
}

Matrix<double> TriangleChiral::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 3;
	tmp(1,0) = 0;
	tmp(0,1) = 0.5;
	tmp(1,1) = sqrt(3.0)/2;
	return tmp;
}

unsigned int TriangleChiral::match_pos_in_ab(Vector<double> const& x) const {
	Vector<double> match(2,0);
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 0; }
	match(0) = 1.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 1; }
	match(0) = 2.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 2; }
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return 3;
}
/*}*/

/*{method needed for checking*/
void TriangleChiral::display_results(){
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
	ps.linked_lines("-",draw_boundary(false),"linecolor=yellow");

	std::complex<double> t;
	unsigned int s0;
	unsigned int s1;
	Matrix<int> nb;
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

			if(my::are_equal(t.imag(),0)){
				arrow = "-";
			} else {
				arrow = "->";
				ps.put((xy0(0)+xy1(0))/2.0,(xy0(1)+xy1(1))/2.0,"\\tiny{"+my::tostring(my::chop(std::arg(-t)/M_PI))+"}");
			}

			ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}
		if(i%3==2){ ps.put(xy0(0)+0.2,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}"); } 
		if(i%3<2){
			unsigned int j(0);
			double flux(0.0);
			do {
				nb = get_neighbourg(s0);
				s1 = nb(i%3+2*j,0);
				flux += std::arg(-H_(s0,s1));
				s0 = s1;
			} while (++j<3);
			flux = my::chop(flux/M_PI);
			if(i%3){ ps.put(xy0(0),xy0(1)+sqrt(3.0)/4.0,"\\tiny{"+my::tostring(flux)+"}"); }
			else   { ps.put((xy0(0)+xy1(0))/2.0,xy0(1)+sqrt(3.0)/4.0,"\\tiny{"+my::tostring(flux)+"}"); }
		}
	}
	ps.end(true,true,true);
}

void TriangleChiral::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="triangle-chiral";
	display_results();

	//compute_H();
	//plot_band_structure();
}
/*}*/
