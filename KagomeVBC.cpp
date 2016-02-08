#include "KagomeVBC.hpp"

KagomeVBC::KagomeVBC(System const& s):
	System(s),
	Kagome<std::complex<double> >(set_ab(),9,"kagome-vbc")
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		system_info_.text("KagomeVBC :");
		system_info_.text(" Each color has the same Hamiltonian.");
		system_info_.text(" 9 sites per unit cell.");
		system_info_.text(" pi-flux through 1/3 of the hexagon and -pi/6-flux through all triangles");
		system_info_.text(" The total flux is null");
	}
}

/*{method needed for running*/
void KagomeVBC::compute_H(){
	H_.set(n_,n_,0);

	double t(-1.0);
	double phi(M_PI/6.0);
	unsigned int s0(0);
	unsigned int s1(0);
	unsigned int ab0(0);
	unsigned int ab1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		ab0 = obs_[0](i,5);
		ab1 = obs_[0](i,6);
		if((ab0==3 && ab1==5) || (ab0==4 && ab1==6) || (ab0==5 && ab1==2) || (ab0==7 && ab1==8) || (ab0==8 && ab1==0)){ H_(s0,s1) = std::polar((obs_[0](i,4)?bc_*t:t),-phi); }
		else { H_(s0,s1) = std::polar((obs_[0](i,4)?bc_*t:t),phi); }
	}
	H_ += H_.conjugate_transpose();
}

void KagomeVBC::create(){
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

Matrix<double> KagomeVBC::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 3.0;
	tmp(1,0) =-sqrt(3.0);
	tmp(0,1) = 0.0;
	tmp(1,1) = 2*sqrt(3.0);
	return tmp;
}

unsigned int KagomeVBC::match_pos_in_ab(Vector<double> const& x) const {
	Vector<double> match(2,0);
	if(my::are_equal(x,match)){ return 0; }
	match(0) = 0.5;
	if(my::are_equal(x,match)){ return 1; }
	match(0) = 0.0;
	match(1) = 0.5;
	if(my::are_equal(x,match)){ return 2; }
	double a(1.0/3.0);
	double b(1.0/6.0);
	match(0) = a;
	match(1) = b; 
	if(my::are_equal(x,match)){ return 3; }
	match(0) += a;
	match(1) += b; 
	if(my::are_equal(x,match)){ return 4; }
	match(0) -= 0.5;
	if(my::are_equal(x,match)){ return 5; }
	match(0) += 2*a;
	match(1) += 2*b;
	if(my::are_equal(x,match)){ return 6; }
	match(0) -= 0.5;
	if(my::are_equal(x,match)){ return 7; }
	match(0) += a;
	match(1) += b;
	if(my::are_equal(x,match)){ return 8; }
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return 9;
}
/*}*/

/*{method needed for checking*/
void KagomeVBC::display_results(){
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
		ps.put(xy0(0)-0.20,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}");
	}
	ps.line("-",boundary_vertex_[0](0),boundary_vertex_[0](1),boundary_vertex_[1](0),boundary_vertex_[1](1),"linecolor=yellow");
	ps.line("-",boundary_vertex_[3](0),boundary_vertex_[3](1),boundary_vertex_[0](0),boundary_vertex_[0](1),"linecolor=yellow");
	ps.end(true,true,true);
}

void KagomeVBC::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="kagome-vbc";
	display_results();

	//plot_band_structure();
}
/*}*/
