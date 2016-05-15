#include "HoneycombChiral.hpp"

HoneycombChiral::HoneycombChiral(System const& s, double const& phi):
	System(s),
	Honeycomb<std::complex<double> >(set_ab(),6,"honeycomb-chiral"),
	phi_(phi)
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		system_info_.text("HoneycombChiral :");
		system_info_.item("Each color has the same Hamiltonian.");
		system_info_.item("There is a flux of "+RST::math(my::tostring(phi)+"\\pi/3") + "per hexagonal plaquette");

		filename_ += "-phi"+my::tostring(phi_);
	 	if(!my::are_equal(phi_,2.0)){ std::cerr<<__PRETTY_FUNCTION__<<" : only well defineed for phi=2 -> phi=2pi/3"<<std::endl; }
	}
}

/*{method needed for running*/
void HoneycombChiral::compute_H(){
	H_.set(n_,n_,0);

	double t(-1.0);
	unsigned int s0(0);
	unsigned int s1(0);
	unsigned int ab0(0);
	unsigned int ab1(0);
	double phi(phi_*M_PI/3.0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		ab0 = obs_[0](i,5);
		ab1 = obs_[0](i,6);
		H_(s0,s1) = (obs_[0](i,4)?bc_*t:t);
		if(ab1==5){
			if(ab0==0){ H_(s0,s1) *= std::polar(1.0,phi); }
			if(ab0==2){ H_(s0,s1) *= std::polar(1.0,-phi); }
		}
	}
	H_ += H_.conjugate_transpose();
}

void HoneycombChiral::create(){
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

void HoneycombChiral::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("phi="+my::tostring(phi_));
		Vector<double> param(1,phi_);

		w.add_header()->title(s,'<');
		w<<param;
		w.add_header()->add(system_info_.get());
	} else { w<<phi_<<" "; }
}

Matrix<double> HoneycombChiral::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 3.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 1.5;
	tmp(1,1) = 1.5*sqrt(3.0);

	//tmp(0,0) = 4.5;
	//tmp(1,0) =-1.5*sqrt(3.0);
	//tmp(0,1) = 0;
	//tmp(1,1) = 2.0*sqrt(3.0);
	return tmp;
}

unsigned int HoneycombChiral::unit_cell_index(Vector<double> const& x) const {
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
	return spuc_;

	//if(my::are_equal(x(1),0.0,eq_prec_,eq_prec_)){
	//if(my::are_equal(x(0),0.0    ,eq_prec_,eq_prec_)){ return 0; }
	//if(my::are_equal(x(0),1.0/3.0,eq_prec_,eq_prec_)){ return 2; }
	//if(my::are_equal(x(0),2.0/3.0,eq_prec_,eq_prec_)){ return 4; }
	//}
	//if(my::are_equal(x(1),0.5,eq_prec_,eq_prec_)){
	//if(my::are_equal(x(0),0.0    ,eq_prec_,eq_prec_)){ return 6; }
	//if(my::are_equal(x(0),1.0/3.0,eq_prec_,eq_prec_)){ return 8; }
	//if(my::are_equal(x(0),2.0/3.0,eq_prec_,eq_prec_)){ return 10; }
	//}
	//if(my::are_equal(x(1),1.0/6.0,eq_prec_,eq_prec_)){
	//if(my::are_equal(x(0),2.0/9.0,eq_prec_,eq_prec_)){ return 1; }
	//if(my::are_equal(x(0),5.0/9.0,eq_prec_,eq_prec_)){ return 3; }
	//if(my::are_equal(x(0),8.0/9.0,eq_prec_,eq_prec_)){ return 5; }
	//}
	//if(my::are_equal(x(1),4.0/6.0,eq_prec_,eq_prec_)){
	//if(my::are_equal(x(0),2.0/9.0,eq_prec_,eq_prec_)){ return 7; }
	//if(my::are_equal(x(0),5.0/9.0,eq_prec_,eq_prec_)){ return 9; }
	//if(my::are_equal(x(0),8.0/9.0,eq_prec_,eq_prec_)){ return 11; }
	//}
	//std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	//return 12;
}
/*}*/

/*{method needed for checking*/
void HoneycombChiral::display_results(){
	compute_H();
	draw_lattice(false,true,ref_(3)?(dir_nn_[4]+dir_nn_[3])*1.5:this->dir_nn_[3]*1.25+this->dir_nn_[4]*0.25);

	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title(RST::math("\\phi")+"="+my::tostring(phi_));
		std::string run_cmd("./mc -s:wf honeycomb-chiral");
		run_cmd += " -u:N " + my::tostring(N_);
		run_cmd += " -u:m " + my::tostring(m_);
		run_cmd += " -u:n " + my::tostring(n_);
		run_cmd += " -i:bc "+ my::tostring(bc_);
		run_cmd += " -d:phi "+ my::tostring(phi_);;
		run_cmd += " -d -u:tmax 10";

		rst_file_->title(title,'-');
		rst_file_->change_text_onclick("run command",run_cmd);

		rst_file_->figure(dir_+filename_+".png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(dir_+filename_+".pdf")+RST::scale("200"));
	}
}

void HoneycombChiral::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="honeycomb-chiral";
	display_results();

	//compute_H();
	//plot_band_structure();
}

void HoneycombChiral::param_fit_therm_limit(std::string& f, std::string& param, std::string& range){
	f="f(x) = a+b*x*x";
	param = "a,b";
	range = "[0:0.015]";
}
/*}*/
