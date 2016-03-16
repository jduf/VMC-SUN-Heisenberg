#include "KagomeFree.hpp"

KagomeFree::KagomeFree(System const& s, Vector<double> const& t):
	System(s),
	Kagome<double>(set_ab(),3,"kagome-free"),
	t_(t)
{
	if(t_.size()==6){
		if(status_==3){ init_lattice(); }
		if(status_==2){
			init_fermionic();

			system_info_.text("KagomeFree :");
			system_info_.item("Each color has the same Hamiltonian.");
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : t must contain 6 values (currently contains "<<t_.size()<<")"<<std::endl; }
}

/*{method needed for running*/
void KagomeFree::compute_H(){
	H_.set(n_,n_,0);

	unsigned int t(0);
	for(unsigned int i(0);i<obs_[0].nlinks(); i++){
		switch(obs_[0](i,5)){
			case 0: { t = (obs_[0](i,3)==1?1:0); } break;
			case 1: { t = (obs_[0](i,3)==2?3:2); } break;
			case 2: { t = (obs_[0](i,3)==1?5:4); } break;
		}
		H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*t_(t):t_(t));
	}
	H_ += H_.transpose();
}

void KagomeFree::create(){
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

void KagomeFree::save_param(IOFiles& w) const {
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
		w.add_header()->add(system_info_.get());
	} else { w<<t_<<" "; }
}

Matrix<double> KagomeFree::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 2.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 1.0;
	tmp(1,1) = sqrt(3.0);
	return tmp;
}

unsigned int KagomeFree::unit_cell_index(Vector<double> const& x) const {
	if(my::are_equal(x(0),0.0,eq_prec_,eq_prec_)){
		if(my::are_equal(x(1),0.0,eq_prec_,eq_prec_)){ return 0; }
		if(my::are_equal(x(1),0.5,eq_prec_,eq_prec_)){ return 2; }
	}
	if(my::are_equal(x(0),0.5,eq_prec_,eq_prec_)){
		if(my::are_equal(x(1),0.0,eq_prec_,eq_prec_)){ return 1; }
	}
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void KagomeFree::lattice(){
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
		if(i%2){ ps.put(xy0(0)-0.20,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}"); }
	}
	ps.end(true,true,true);
}

void KagomeFree::display_results(){
	lattice();

	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title("t=(");
		std::string run_cmd("./mc -s:wf kagome-free");
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

void KagomeFree::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="kagome-free";
	display_results();

	//compute_H();
	//plot_band_structure();
}
/*}*/
