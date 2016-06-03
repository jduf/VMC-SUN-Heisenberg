#include "TriangleT3x2.hpp"

TriangleT3x2::TriangleT3x2(System const& s, Vector<double> const& t):
	System(s),
	Triangle<double>(set_ab(),6,"triangle-T3x2"),
	t_(t)
{
	if(t_.size()==18){
		if(status_==3){ init_lattice(); }
		if(status_==2){
			init_fermionic();

			system_info_.text("TriangleT3x2 :");
			system_info_.item("Each color has a the same Hamiltonian.");

			filename_ += "-t";
			for(unsigned int i(0);i<t_.size();i++){
				filename_ += ((t_(i)>=0)?"+":"")+my::tostring(t_(i));
			}
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : t must contain 18 values (currently contains "<<t_.size()<<")"<<std::endl; }
}

/*{method needed for running*/
void TriangleT3x2::compute_H(){
	H_.set(n_,n_,0);

	unsigned int s0(0);
	unsigned int s1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);

		H_(s0,s1) = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)+3*obs_[0](i,5));
	}
	H_ += H_.transpose();
}

void TriangleT3x2::create(){
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

void TriangleT3x2::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("t=(");

		for(unsigned int i(0);i<t_.size()-1;i++){
			s += my::tostring(t_(i))+",";
		}
		s += my::tostring(t_.back())+") ";

		w.add_header()->title(s,'<');
		w<<t_;
		w.add_header()->add(system_info_.get());
	} else { w<<t_<<" "; }
}

Matrix<double> TriangleT3x2::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 3.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 1.0;
	tmp(1,1) = sqrt(3.0);
	return tmp;
}

unsigned int TriangleT3x2::unit_cell_index(Vector<double> const& x) const {
	if(my::are_equal(x(1),0.0,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.0    ,eq_prec_,eq_prec_)){ return 0; }
		if(my::are_equal(x(0),1.0/3.0,eq_prec_,eq_prec_)){ return 1; }
		if(my::are_equal(x(0),2.0/3.0,eq_prec_,eq_prec_)){ return 2; }
	}
	if(my::are_equal(x(1),0.5,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.0    ,eq_prec_,eq_prec_)){ return 3; }
		if(my::are_equal(x(0),1.0/3.0,eq_prec_,eq_prec_)){ return 4; }
		if(my::are_equal(x(0),2.0/3.0,eq_prec_,eq_prec_)){ return 5; }
	}
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void TriangleT3x2::display_results(){
	compute_H();
	draw_lattice(false,true,ref_(3)?(dir_nn_[4]+dir_nn_[3])*1.5:dir_nn_[3]*1.75+dir_nn_[4]*0.25);

	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title("t=(");
		std::string run_cmd("./mc -s:wf triangle-T3x2");
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

void TriangleT3x2::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="triangle-T3x2";
	display_results();

	//compute_H();
	//plot_band_structure();
}
/*}*/
