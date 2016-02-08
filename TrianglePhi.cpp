#include "TrianglePhi.hpp"

TrianglePhi::TrianglePhi(System const& s, double const& phi):
	System(s),
	Triangle<std::complex<double> >(set_ab(),1,"triangle-phi"),
	phi_(phi)
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		system_info_.text("SquarePhiflux :");
		system_info_.text(" Each color has a different Hamiltonian.");
		system_info_.text(" Each neighbouring triangle has a flux of opposite sign.");

		filename_ += "-phi" + my::tostring(phi_);
	}
}

/*{method needed for running*/
void TrianglePhi::compute_H(){
	H_.set(n_,n_,0);

	double unit_flux(2*M_PI*m_/N_);
	double t(-1.0);
	unsigned int s0(0);
	unsigned int s1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		H_(s0,s1) = (obs_[0](i,4)?bc_*t:t);
		if(obs_[0](i,3)==1){ H_(s0,s1) *= std::polar(1.0,phi_*unit_flux); }
	}
	H_ += H_.conjugate_transpose();
}

void TrianglePhi::create(){
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

void TrianglePhi::save_param(IOFiles& w) const {
	std::string s("phi=("+my::tostring(phi_)+")");
	Vector<double> param(1,phi_);

	w.add_header()->title(s,'<');
	w<<param;
	GenericSystem<std::complex<double> >::save_param(w);
}

Matrix<double> TrianglePhi::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 1;
	tmp(1,0) = 0;
	tmp(0,1) = 0.5;
	tmp(1,1) = sqrt(3.0)/2.0;
	return tmp;
}
/*}*/

/*{method needed for checking*/
void TrianglePhi::display_results(){
	compute_H();

	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1mm");
	std::string arrow("-");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(info_+path_+dir_,filename_);
	ps.begin(-20,-20,20,20,filename_);
	ps.polygon(cluster_vertex_,"linecolor=green");
	ps.polygon(draw_unit_cell(),"linecolor=black");

	std::complex<double> t;
	double unit_flux(2*M_PI*m_/N_);
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

			arrow = "-";
			if(t.imag()>0){ arrow = "-"+std::string(std::arg(t)/unit_flux,'>'); }
			if(t.imag()<0){ arrow = std::string(-std::arg(t)/unit_flux,'<')+"-"; }

			linewidth = my::tostring(std::abs(t))+"mm";
			ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}
		if(i%3==2){ ps.put(xy0(0)+0.2,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}"); }
	}
	ps.line("-",boundary_vertex_[0](0),boundary_vertex_[0](1),boundary_vertex_[1](0),boundary_vertex_[1](1),"linecolor=yellow");
	ps.line("-",boundary_vertex_[3](0),boundary_vertex_[3](1),boundary_vertex_[0](0),boundary_vertex_[0](1),"linecolor=yellow");
	ps.end(true,true,true);
}

void TrianglePhi::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="triangle-phi";
	display_results();
}
/*}*/
