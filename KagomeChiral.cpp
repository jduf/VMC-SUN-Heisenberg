#include "KagomeChiral.hpp"

KagomeChiral::KagomeChiral(System const& s, double const& phi):
	System(s),
	Kagome<std::complex<double> >(set_ab(),6,"kagome-chiral"),
	phi_(phi)
{
	if(phi<=3.0){
		if(status_==3){ init_lattice(); }
		if(status_==2){
			init_fermionic();

			system_info_.text("KagomeChiral :");
			system_info_.item("Each color has the same Hamiltonian.");
			system_info_.item("6 sites per unit cell.");
			system_info_.item("Flux of "+RST::math(my::tostring(phi)+"\\pi/3")+" per plaquette.");
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : the flux per plaquette shouldn't be bigger than pi (phi<=3)"<<std::endl; }
}

/*{method needed for running*/
void KagomeChiral::compute_H(){
	H_.set(n_,n_,0);

	double t(-1.0);
	double phi(phi_*M_PI/3.0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		switch(obs_[0](i,5)){
			case 0:{ H_(obs_[0](i,0),obs_[0](i,1)) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)==0?-phi:0.0)); }break;
			case 1:{ H_(obs_[0](i,0),obs_[0](i,1)) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)==2?phi:0.0)); }break;
			case 2:{ H_(obs_[0](i,0),obs_[0](i,1)) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)==0?3.0*phi:0.0)); }break;
			case 3:{ H_(obs_[0](i,0),obs_[0](i,1)) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)==0?-phi:0.0)); }break;
			case 4:{ H_(obs_[0](i,0),obs_[0](i,1)) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)==2?-2.0*phi:0.0)); }break;
			case 5:{ H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*t:t); }break;
		}
	}
	H_ += H_.conjugate_transpose();
}

void KagomeChiral::create(){
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

void KagomeChiral::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("phi="+my::tostring(phi_));
		Vector<double> param(1,phi_);

		w.add_header()->title(s,'<');
		w<<param;
		w.add_header()->add(system_info_.get());
	} else { w<<phi_<<" "; }
}

Matrix<double> KagomeChiral::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 2.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 0.0;
	tmp(1,1) = 2.0*sqrt(3.0);
	return tmp;
}

unsigned int KagomeChiral::unit_cell_index(Vector<double> const& x) const {
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
void KagomeChiral::display_results(){
	compute_H();
	draw_lattice();

	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title(RST::math("\\phi")+"="+my::tostring(phi_));
		std::string run_cmd("./mc -s:wf kagome-chiral");
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

void KagomeChiral::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="kagome-chiral";
	display_results();

	//compute_H();
	//plot_band_structure();
}
/*}*/
