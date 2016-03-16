#include "KagomePlaquette6B.hpp"

KagomePlaquette6B::KagomePlaquette6B(System const& s, double const& td):
	System(s),
	Kagome<double>(set_ab(),6,"kagome-plaquette6B"),
	td_(td)
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		filename_ += ((td_>=0)?"-td+":"-td") + my::tostring(td_);

		system_info_.text("KagomePlaquette6B :");
		system_info_.item("Each color has the same Hamiltonian.");
		system_info_.item("6 sites per unit cell.");
		if(td_<0.0){ system_info_.item(RST::math("(\\pi;\\pi,\\pi)")+"-flux"); }
		if(td_>0.0){ system_info_.item(RST::math("(0;0,\\pi)")+"-flux"); }
	}
}

/*{method needed for running*/
void KagomePlaquette6B::compute_H(){
	H_.set(n_,n_,0);

	double th(-1.0);
	double t(0.0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		t=0.0;
		switch(obs_[0](i,5)){
			case 0: { t = th; } break;
			case 1: { t = td_; } break;
			case 2: { t = (obs_[0](i,3)?th:td_); } break;
			case 3: { t = (obs_[0](i,3)?th:-th);  } break;
			case 4: { t = td_; } break;
			case 5: { t = (obs_[0](i,3)?-th:td_); } break;
		}
		H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*t:t);
	}
	H_ += H_.transpose();
}

void KagomePlaquette6B::create(){
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

void KagomePlaquette6B::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("td="+my::tostring(td_));
		Vector<double> param(1,td_);

		w.add_header()->title(s,'<');
		w<<param;
		w.add_header()->add(system_info_.get());
	} else { w<<td_<<" "; }
}

Matrix<double> KagomePlaquette6B::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 2.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 0.0;
	tmp(1,1) = 2.0*sqrt(3.0);
	return tmp;
}

unsigned int KagomePlaquette6B::unit_cell_index(Vector<double> const& x) const {
	if(my::are_equal(x(1),0.0,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.0,eq_prec_,eq_prec_)){ return 0; }
		if(my::are_equal(x(0),0.5,eq_prec_,eq_prec_)){ return 2; }
	}
	if(my::are_equal(x(1),0.25,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.25,eq_prec_,eq_prec_)){ return 1; }
	}
	if(my::are_equal(x(1),0.5,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.0,eq_prec_,eq_prec_)){ return 5; }
		if(my::are_equal(x(0),0.5,eq_prec_,eq_prec_)){ return 3; }
	}
	if(my::are_equal(x(1),0.75,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.75,eq_prec_,eq_prec_)){ return 4; }
	}
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void KagomePlaquette6B::lattice(){
	Vector<unsigned int> o(3,0);
	for(unsigned int i(1);i<obs_.size();i++){
		switch(obs_[i].get_type()){
			case 1:{ o(0)=i; }break;//bond energy
			case 2:{ o(1)=i; }break;//long range correlation
			case 3:{ o(2)=i; }break;//color occupation
		}
	}
	compute_H();

	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(info_+path_+dir_,filename_);
	ps.begin(-20,-20,20,20,filename_);
	ps.polygon(cluster_vertex_,"linecolor=green");
	Matrix<double> uc(draw_unit_cell(-1.0,-sqrt(3.0)/4.0));
	ps.polygon(uc,"linecolor=black");
	ps.linked_lines("-",draw_boundary(false),"linecolor=yellow");

	if(o(1)){ draw_long_range_correlation(ps,obs_[o(1)]); }
	double t;
	unsigned int s0;
	unsigned int s1;
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		xy0 = x_[s0];

		s1 = obs_[0](i,1);
		xy1 = x_[s1];

		t = H_(s0,s1);
		linewidth = my::tostring(t)+"mm";
		if(o(0) || o(2)){
			if(my::in_polygon(uc.row(),uc.ptr(),uc.ptr()+uc.row(),xy0(0),xy0(1))){
				if(o(0)){ t = obs_[o(0)][obs_[0](i,2)].get_x(); }
				if(i%2 && o(2)){
					Vector<double> p(N_);
					for(unsigned int j(0);j<N_;j++){ p(j) = obs_[o(2)][j+N_*obs_[0](i,5)].get_x(); }
					ps.pie(xy0(0),xy0(1),p,0.2,"chartColor=color");
				}
				linewidth = my::tostring(std::abs(t))+"mm";
			} else if(my::in_polygon(uc.row(),uc.ptr(),uc.ptr()+uc.row(),xy1(0),xy1(1))){ t = 0; }
		}
		if(std::abs(t)>1e-4){
			if((xy0-xy1).norm_squared()>1.0001){
				linestyle = "dashed";
				xy1 = (xy0+dir_nn_[obs_[0](i,3)]).chop();
				ps.put(xy1(0)-0.20,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
			} else{ linestyle = "solid"; }

			if(t>0){ color = "blue"; }
			else   { color = "red"; }
			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}
		if(i%2){
			ps.put(xy0(0)-0.20,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}"); 
			ps.put(xy0(0)-0.20,xy0(1)-0.15,"\\textcolor{green}{\\tiny{"+my::tostring(obs_[0](i,5))+"}}"); 
		}
	}
	ps.end(true,true,true);
}

void KagomePlaquette6B::display_results(){
	lattice();

	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title("t="+my::tostring(td_));
		std::string run_cmd("./mc -s:wf kagome-plaquette6B");
		run_cmd += " -u:N " + my::tostring(N_);
		run_cmd += " -u:m " + my::tostring(m_);
		run_cmd += " -u:n " + my::tostring(n_);
		run_cmd += " -i:bc "+ my::tostring(bc_);
		run_cmd += " -d:td "+ my::tostring(td_);
		run_cmd += " -d -u:tmax 10";

		rst_file_->title(title,'-');
		rst_file_->change_text_onclick("run command",run_cmd);

		rst_file_->figure(dir_+filename_+".png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(dir_+filename_+".pdf")+RST::scale("200"));
	}
}

void KagomePlaquette6B::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="kagome-plaquette6B";
	display_results();
}
/*}*/
