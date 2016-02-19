#include "Honeycomb0pp.hpp"

Honeycomb0pp::Honeycomb0pp(System const& s, double const& td, unsigned int const& fc):
	System(s),
	Honeycomb<double>(set_ab(),6,"honeycomb-0pp"),
	td_(td),
	fc_(fc)
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		system_info_.text("Honeycomb0pp :");
		system_info_.item("Each color has the same Hamiltonian.");
		system_info_.item("6 sites per unit cell");
		system_info_.item("1/3 of the hexagons with 0-flux.");
		system_info_.item("If td<0, each 0-flux hexagon is surrounded by pi-flux hexagons, 0-flux otherwise.");
		system_info_.item("th is set to -1");

		filename_ += "-td" + my::tostring(td_)+"-fc" + my::tostring(fc_);
	}
}

/*{method needed for running*/
void Honeycomb0pp::compute_H(){
	H_.set(n_,n_,0);

	double th(-1.0);
	unsigned int s0(0);
	unsigned int s1(0);
	unsigned int ab0(0);
	unsigned int ab1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		ab0 = obs_[0](i,5);
		ab1 = obs_[0](i,6);
		switch(fc_){
			case 0:
				{
					if((ab0==0 && ab1==1) || (ab0==2 && ab1==3) || (ab0==4 && ab1==5)){ H_(s0,s1) = (obs_[0](i,4)?bc_*td_:td_); } 
					else { H_(s0,s1) = (obs_[0](i,4)?bc_*th:th); }
				}break;
			case 1:
				{
					if((ab0==4 && ab1==3) || (ab0==2 && ab1==1) || (ab0==0 && ab1==5)){ H_(s0,s1) = (obs_[0](i,4)?bc_*td_:td_); } 
					else { H_(s0,s1) = (obs_[0](i,4)?bc_*th:th); }
				}break;
			case 2:
				{
					if((ab0==2 && ab1==5) || (ab0==0 && ab1==3) || (ab0==4 && ab1==1)){ H_(s0,s1) = (obs_[0](i,4)?bc_*td_:td_); } 
					else { H_(s0,s1) = (obs_[0](i,4)?bc_*th:th); }
				}break;
			default:
				{ std::cerr<<__PRETTY_FUNCTION__<<" unknown flux configuration 3>fc="<<fc_<<std::endl; }
		}
	}
	H_ += H_.transpose();
}

void Honeycomb0pp::create(){
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

void Honeycomb0pp::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("td="+my::tostring(td_)+", th=-1, fc="+my::tostring(fc_));
		Vector<double> param(2);
		param(0) = td_;
		param(1) = fc_;

		w.add_header()->title(s,'<');
		w<<param;
		GenericSystem<double>::save_param(w);
	} else { w<<td_<<" "<<fc_<<" "; }
}

Matrix<double> Honeycomb0pp::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 3.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 1.5;
	tmp(1,1) = 1.5*sqrt(3.0);
	return tmp;
}

unsigned int Honeycomb0pp::match_pos_in_ab(Vector<double> const& x) const {
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
void Honeycomb0pp::lattice(){
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
	ps.linked_lines("-",draw_boundary(false),"linecolor=yellow");

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
			xy0 -=  dir_nn_[obs_[0](i,3)]*0.2;
			xy1 -=  dir_nn_[obs_[0](i,3)]*0.2;
			ps.put(xy0(0),xy0(1),"\\tiny{"+my::tostring(s0)+"}"); 
			xy0 -=  dir_nn_[obs_[0](i,3)]*0.2;
			xy1 -=  dir_nn_[obs_[0](i,3)]*0.2;
			ps.put(xy0(0),xy0(1),"\\textcolor{green}{\\tiny{"+my::tostring(obs_[0](i,5))+"}}"); 

			xy0 +=  dir_nn_[obs_[0](i,3)]*0.6;
			xy1 +=  dir_nn_[obs_[0](i,3)]*0.6;
			ps.put(xy1(0),xy1(1),"\\tiny{"+my::tostring(s1)+"}");
			xy0 +=  dir_nn_[obs_[0](i,3)]*0.2;
			xy1 +=  dir_nn_[obs_[0](i,3)]*0.2;
			ps.put(xy1(0),xy1(1),"\\textcolor{green}{\\tiny{"+my::tostring(obs_[0](i,6))+"}}");
		}
	}
	ps.end(true,true,true);
}

void Honeycomb0pp::display_results(){
	lattice();
	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string run_cmd("./mc -s:wf ladder-free");
		run_cmd += " -u:N " + my::tostring(N_);
		run_cmd += " -u:m " + my::tostring(m_);
		run_cmd += " -u:n " + my::tostring(n_);
		run_cmd += " -i:bc "+ my::tostring(bc_);
		run_cmd += " -d:td " + my::tostring(td_); 
		run_cmd += " -u:fc " + my::tostring(fc_); 
		run_cmd += " -d:Jp 1 -u:tmax 10 -d";
		rst_file_->change_text_onclick("run command",run_cmd);

		rst_file_->figure(dir_+filename_+"-pstricks.png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(dir_+filename_+"-pstricks.pdf")+RST::scale("200"));
	}
}

void Honeycomb0pp::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="honeycomb-0pp";
	display_results();

	//plot_band_structure();
}
/*}*/
