#include "SquareDimerizedBar.hpp"

SquareDimerizedBar::SquareDimerizedBar(System const& s, Vector<double> const& t):
	System(s),
	Square<double>(set_ab(),4,"square-dimerizedbar"),
	t_(t)
{
	if(t.size()==8){
		if(status_==3){ init_lattice(); }
		if(status_==2){
			init_fermionic();

			system_info_.text("SquareDimerizedBar :");
			system_info_.item("Each color has the same Hamiltonian.");
			system_info_.item("4 sites in a 2x2 unit cell");
			system_info_.item("The dimerization occurs on bonds stacked in columns with strongest "+RST::math("\\vert t\\vert"));
			system_info_.item(RST::math("1=t_0\\geq t_i\\geq 0")+" (signs are set in H)");

			filename_ += "-t";
			for(unsigned int i(0);i<t_.size();i++){
				filename_ += ((t_(i)>=0)?"+":"")+my::tostring(t_(i));
			}
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : t must contain 8 values"<<std::endl; }
}

/*{method needed for running*/
void SquareDimerizedBar::compute_H(){
	H_.set(n_,n_,0);

	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,5)==1?-1:1)*(obs_[0](i,4)?bc_:1)*t_(2*obs_[0](i,5)+obs_[0](i,3));
	}
	H_ += H_.transpose();
}

void SquareDimerizedBar::create(){
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

void SquareDimerizedBar::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("t=(");
		Vector<double> param(t_.size());

		for(unsigned int i(0);i<t_.size()-1;i++){
			param(i) = t_(i);
			s += my::tostring(t_(i))+",";
		}
		param(t_.size()-1) = t_.back();
		s += my::tostring(t_.back())+")";

		w.add_header()->title(s,'<');
		w<<param;
		GenericSystem<double>::save_param(w);
	} else { w<<t_<<" "; }
}

Matrix<double> SquareDimerizedBar::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 2.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 0.0;
	tmp(1,1) = 2.0;
	return tmp;
}

unsigned int SquareDimerizedBar::match_pos_in_ab(Vector<double> const& x) const {
	unsigned int i(0);
	if(my::are_equal(x(0),0.5,eq_prec_,eq_prec_)){ i+=1; }
	if(my::are_equal(x(1),0.5,eq_prec_,eq_prec_)){ i+=2; }
	return i;
}
/*}*/

/*{method needed for checking*/
void SquareDimerizedBar::lattice(){
	compute_H();

	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1mm");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(info_+path_+dir_,filename_);
	ps.begin(-20,-20,20,20,filename_);
	ps.polygon(cluster_vertex_,"linecolor=green");
	Matrix<double> uc(draw_unit_cell(0.5,0.5));
	ps.polygon(uc,"linecolor=black");
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
		if(obs_.size()>1){
			if(my::in_polygon(uc.row(),uc.ptr(),uc.ptr()+uc.row(),xy0(0),xy0(1))){ t = obs_[1][i%4].get_x(); }
			else if(my::in_polygon(uc.row(),uc.ptr(),uc.ptr()+uc.row(),xy1(0),xy1(1))){ t = 0; }
		}
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

		if(i%2){ ps.put(xy0(0)+0.2,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}"); }
	}
	ps.end(true,true,true);
}

void SquareDimerizedBar::display_results(){
	lattice();

	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title("t=(");
		std::string run_cmd("./mc -s:wf square-dimerizedbar");
		run_cmd += " -u:N " + my::tostring(N_);
		run_cmd += " -u:m " + my::tostring(m_);
		run_cmd += " -u:n " + my::tostring(n_);
		run_cmd += " -i:bc "+ my::tostring(bc_);
		run_cmd += " -d:t ";
		for(unsigned int i(0);i<t_.size()-1;i++){
			title   += my::tostring(t_(i)) + ","; 
			run_cmd += my::tostring(t_(i)) + ","; 
		}
		title   += my::tostring(t_.back()) + ")";
		run_cmd += my::tostring(t_.back()) + " -d -u:tmax 10";

		rst_file_->title(title,'-');
		rst_file_->change_text_onclick("run command",run_cmd);

		rst_file_->figure(dir_+filename_+".png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(dir_+filename_+".pdf")+RST::scale("200"));
	}
}

void SquareDimerizedBar::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="square-dimerizedbar";
	display_results();

	//plot_band_structure();
}
/*}*/
