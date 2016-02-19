#include "TriangleMu.hpp"

TriangleMu::TriangleMu(System const& s, double const& mu):
	System(s),
	Triangle<double>(set_ab(),3,"triangle-mu"),
	mu_(mu)
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();
		same_wf_ = false;

		system_info_.text("TriangleMu :");
		system_info_.item("Each color has a different Hamiltonian.");

		filename_ += "-mu"+std::string(mu_>=0?"+":"")+my::tostring(mu_);
	}
}

/*{method needed for running*/
void TriangleMu::compute_H(unsigned int const& c){
	H_.set(n_,n_,0);

	double t(-1.0);
	unsigned int s0(0);
	unsigned int s1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		H_(s0,s1) = (obs_[0](i,4)?bc_*t:t);
		if((unsigned int)(obs_[0](i,5))==c%spuc_){ H_(s0,s0) = mu_/2; }
	}
	H_ += H_.transpose();
}

void TriangleMu::create(){
	for(unsigned int c(0);c<N_;c++){
		status_ = 2;
		compute_H(c);
		diagonalize(true);
		if(status_==1){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = H_(i,j);
				}
			}
		}
	}
}

void TriangleMu::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("mu=("+my::tostring(mu_)+")");
		Vector<double> param(1,mu_);

		w.add_header()->title(s,'<');
		w<<param;
		GenericSystem<double>::save_param(w);
	} else { w<<mu_<<" "; }
}

Matrix<double> TriangleMu::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 1.5;
	tmp(1,0) =-sqrt(3.0)/2;
	tmp(0,1) = 1.5;
	tmp(1,1) = sqrt(3.0)/2;
	return tmp;
}

unsigned int TriangleMu::unit_cell_index(Vector<double> const& x) const {
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
/*}*/

/*{method needed for checking*/
void TriangleMu::display_results(){
	compute_H(0);

	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(info_+path_+dir_,filename_);
	ps.begin(-20,-20,20,20,filename_);
	ps.polygon(cluster_vertex_,"linecolor=green");
	ps.polygon(draw_unit_cell(),"linecolor=black");
	ps.linked_lines("-",draw_boundary(false),"linecolor=yellow");

	double t;
	double mu;
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

			if(t>0){ color = "blue"; }
			else   { color = "red"; }
			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);

			mu = H_(s0,s0);
			if(std::abs(mu)>1e-4){
				if(mu>0){ color = "cyan"; }
				else    { color = "magenta"; }
				ps.circle(xy0,std::abs(mu),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
			}
		}
		if(i%3==2){ ps.put(xy0(0)+0.2,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}"); }
	}
	ps.end(true,true,true);
}

void TriangleMu::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="triangle-mu";
	display_results();

	//plot_band_structure();
}
/*}*/
