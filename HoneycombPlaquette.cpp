#include "HoneycombPlaquette.hpp"

HoneycombPlaquette::HoneycombPlaquette(System const& s, Vector<double> const& t):
	System(s),
	Honeycomb<double>(set_ab(),6,"honeycomb-plaquette"),
	t_(t)
{
	if(status_==2){
		init_lattice();
		init_fermionic();

		filename_ += "-t";
		for(unsigned int i(0);i<t_.size();i++){
			filename_ += ((t_(i)>=0)?"+":"")+my::tostring(t_(i));
		}
	}
}

/*{method needed for running*/
void HoneycombPlaquette::compute_H(){
	H_.set(n_,n_,0);

	unsigned int s0(0);
	unsigned int s1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		if( std::abs(obs_[0](i,5)-obs_[0](i,6)) == 3 ){ H_(s0,s1) = (obs_[0](i,4)?bc_*t_(1):t_(1)); }
		else { H_(s0,s1) = (obs_[0](i,4)?bc_*t_(0):t_(0)); }
	}
	H_ += H_.transpose();
}

void HoneycombPlaquette::create(){
	compute_H();
	diagonalize(true);
	for(unsigned int c(0);c<N_;c++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_(c);j++){
				EVec_[c](i,j) = H_(i,j);
			}
		}
	}
}

void HoneycombPlaquette::save_param(IOFiles& w) const {
	std::string s("t=(");

	for(unsigned int i(0);i<t_.size()-1;i++){ 
		s += my::tostring(t_(i))+",";
	}
	s += my::tostring(t_.back())+")";

	w.add_header()->title(s,'<');
	w<<t_;
	GenericSystem<double>::save_param(w);
}

Matrix<double> HoneycombPlaquette::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 3.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 1.5;
	tmp(1,1) = 1.5*sqrt(3.0);
	return tmp;
}

unsigned int HoneycombPlaquette::match_pos_in_ab(Vector<double> const& x) const {
	Vector<double> match(2,0);
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 0; }
	match(0) = 1.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 1; }
	match(1) = 1.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 2; }
	match(0) = 0;
	match(1) = 2.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 3; }
	match(0) = 2.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 4; }
	match(1) = 1.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 5; }
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return 6;
}
/*}*/

/*{method needed for checking*/
void HoneycombPlaquette::lattice(){
	compute_H();

	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(info_+path_+dir_,filename_);
	ps.begin(-20,-20,20,20,filename_);
	ps.polygon(cluster_vertex_,"linecolor=green");
	ps.polygon(draw_unit_cell(),"linecolor=black");

	double t;
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
				if(obs_[0](i,3)){ ps.put(xy1(0)+0.2,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}"); }
			} else { linestyle = "solid"; }

			if(t>0){ color = "blue";}
			else   { color = "red"; }
			linewidth = my::tostring(std::abs(t))+"mm";
			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}
		if(!(i%3)){ 
			ps.put(xy1(0)+0.10,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
			ps.put(xy0(0)+0.10,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}"); 
		}
	}
	ps.line("-",boundary_vertex_[0](0),boundary_vertex_[0](1),boundary_vertex_[1](0),boundary_vertex_[1](1),"linecolor=yellow");
	ps.line("-",boundary_vertex_[3](0),boundary_vertex_[3](1),boundary_vertex_[0](0),boundary_vertex_[0](1),"linecolor=yellow");
	ps.end(true,true,true);
}

void HoneycombPlaquette::display_results(){
	lattice();
	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title("t=(");
		std::string run_cmd("./mc -s:wf ladder-free");
		run_cmd += " -u:N " + my::tostring(N_);
		run_cmd += " -u:m " + my::tostring(m_);
		run_cmd += " -u:n " + my::tostring(n_);
		run_cmd += " -i:bc "+ my::tostring(bc_);
		run_cmd += " -d:t ";
		for(unsigned int i(0);i<t_.size()-1;i++){
			title   += my::tostring(t_(i)) + ","; 
			run_cmd += my::tostring(t_(i)) + ","; 
		}
		title   += my::tostring(t_.back()) + "),";
		run_cmd += " -d:Jp 1 -u:tmax 10 -d";

		rst_file_->title(title,'-'); 
		rst_file_->change_text_onclick("run command",run_cmd);

		rst_file_->figure(dir_+filename_+"-pstricks.png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(dir_+filename_+"-pstricks.pdf")+RST::scale("200"));
	}
}

void HoneycombPlaquette::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="honeycomb-plaquette";
	display_results();
}
/*}*/
