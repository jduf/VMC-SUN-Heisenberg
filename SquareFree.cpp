#include "SquareFree.hpp"

SquareFree::SquareFree(System const& s, Vector<double> const& t, Vector<double> const& mu):
	System(s),
	Square<double>(set_ab(),5,"square-free"),
	t_(t),
	mu_(mu)
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		bool need_compute_additional_links(true);
		for(unsigned int i(0);i<obs_.size();i++){
			if(obs_[i].get_type() == 4){ need_compute_additional_links = false; i=obs_.size(); }
		}
		if(need_compute_additional_links){ init_additional_links(); }
		same_wf_ = false;

		system_info_.text("SquareFree :");
		system_info_.text(" Each colors has a different Hamiltonian.");
		system_info_.text(" There is an additional second neighbour hopping for every 1/k sites.");

		filename_ += "-t";
		for(unsigned int i(0);i<t_.size();i++){
			filename_ += ((t_(i)>=0)?"+":"")+my::tostring(t_(i));
		}
		filename_ += "-mu";
		for(unsigned int i(0);i<mu_.size();i++){
			filename_ += ((mu_(i)>=0)?"+":"")+my::tostring(mu_(i));
		}
	}
}

/*{method needed for running*/
void SquareFree::init_additional_links(){
	Vector<double> x;
	Matrix<int> tmp(n_*2,3);
	for(unsigned int i(0);i<n_;i++){
		x = x_[i]+dir_nn_[0]+dir_nn_[1]*2.0;
		tmp(2*i,2) = handle_boundary(x_[i],x);
		tmp(2*i,0) = i;
		tmp(2*i,1) = find_index(x);

		x = x_[i]-dir_nn_[1]+dir_nn_[0]*2.0;
		tmp(2*i+1,2) = handle_boundary(x_[i],x);
		tmp(2*i+1,0) = i;
		tmp(2*i+1,1) = find_index(x);
	}
	obs_.push_back(Observable("Additional links",4,0,tmp));
}

void SquareFree::compute_H(unsigned int const& c){
	H_.set(n_,n_,0);

	unsigned int s0(0);
	unsigned int s1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		H_(s0,s1) = (obs_[0](i,4)?bc_*t_(0):t_(0));
		if((unsigned int)(obs_[0](i,5))==c%spuc_){
			if(obs_[0](i,3)){ H_(s0,s0) = mu_(0)/2; }
			H_(obs_[1](i,0),obs_[1](i,1)) = (obs_[1](i,2)?bc_*t_(1):t_(1));
		}
	}
	H_ += H_.transpose();
}

void SquareFree::create(){
	for(unsigned int c(0);c<N_;c++){
		compute_H(c);
		status_=2;
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

void SquareFree::save_param(IOFiles& w) const {
	std::string s("t=(");
	Vector<double> param(t_.size()+mu_.size());

	for(unsigned int i(0);i<t_.size()-1;i++){
		param(i) = t_(i);
		s += my::tostring(t_(i))+",";
	}
	param(t_.size()-1) = t_.back();
	s += my::tostring(t_.back())+") "+RST::math("\\mu")+"=(";

	for(unsigned int i(0);i<mu_.size()-1;i++){
		param(i+t_.size()) = mu_(i);
		s += my::tostring(mu_(i))+",";
	}
	param.back() = mu_.back();
	s += my::tostring(mu_.back())+")";

	w.add_header()->title(s,'<');
	w<<param;
	GenericSystem<double>::save_param(w);
}

Matrix<double> SquareFree::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 2.0;
	tmp(1,0) =-1.0;
	tmp(0,1) = 1.0;
	tmp(1,1) = 2.0;
	return tmp;
}

unsigned int SquareFree::match_pos_in_ab(Vector<double> const& x) const {
	Vector<double> match(2,0);
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
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return 5;
}
/*}*/

/*{method needed for checking*/
void SquareFree::display_results(){
	compute_H(0);
	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1mm");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(info_+path_+dir_,filename_);
	ps.begin(-20,-20,20,20,filename_);
	ps.polygon(cluster_vertex_,"linecolor=green");
	ps.polygon(draw_unit_cell(),"linecolor=black");

	double t;
	double mu;
	unsigned int s0;
	unsigned int s1;
	for(unsigned int i(0);i<obs_[1].nlinks();i++){
		s0 = obs_[1](i,0);
		xy0 = x_[s0];

		s1 = obs_[1](i,1);
		xy1 = x_[s1];

		t = H_(s0,s1);
		if(std::abs(t)>1e-4){
			if((xy0-xy1).norm_squared()>5.001){
				linestyle = "dashed";
				if(i%2){ xy1=xy0-dir_nn_[1]+dir_nn_[0]*2.0; }
				else   { xy1=xy0+dir_nn_[0]+dir_nn_[1]*2.0; }
				ps.put(xy1(0)+0.2,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
			} else { linestyle = "solid"; }

			if(t>0){ color = "blue"; }
			else   { color = "red"; }
			linewidth=my::tostring(std::abs(t))+"mm";
			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}
	}
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
			linewidth=my::tostring(std::abs(t))+"mm";
			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}

		mu = H_(s0,s0);
		if(std::abs(mu)>1e-4){
			if(mu<0){ color = "green"; }
			else    { color = "cyan"; }
			ps.circle(xy0,std::abs(mu),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
		}

		if(i%2){ ps.put(xy0(0)+0.2,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}"); }
	}
	ps.line("-",boundary_vertex_[0](0),boundary_vertex_[0](1),boundary_vertex_[1](0),boundary_vertex_[1](1),"linecolor=yellow");
	ps.line("-",boundary_vertex_[3](0),boundary_vertex_[3](1),boundary_vertex_[0](0),boundary_vertex_[0](1),"linecolor=yellow");
	ps.end(true,true,true);
}

void SquareFree::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="square-free";
	display_results();

	//plot_band_structure();
}
/*}*/
