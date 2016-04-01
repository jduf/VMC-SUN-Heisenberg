#include "TrianglePhi.hpp"

TrianglePhi::TrianglePhi(System const& s, double const& phi):
	System(s),
	Triangle<std::complex<double> >(set_ab(),1,"triangle-phi"),
	phi_(phi)
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		system_info_.text("SquarePhiflux :");
		system_info_.item("Each color has the same Hamiltonian.");
		system_info_.item("Each neighbouring triangle has a flux of opposite sign.");

		filename_ += "-phi" + my::tostring(phi_);
	}
}

/*{method needed for running*/
void TrianglePhi::compute_H(){
	H_.set(n_,n_,0);

	double unit_flux(2*M_PI*m_/N_);
	double t(-1.0);
	unsigned int s0(0);
	unsigned int s1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		H_(s0,s1) = (obs_[0](i,4)?bc_*t:t);
		if(obs_[0](i,3)==1){ H_(s0,s1) *= std::polar(1.0,phi_*unit_flux); }
	}
	H_ += H_.conjugate_transpose();
}

void TrianglePhi::create(){
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

void TrianglePhi::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("phi=("+my::tostring(phi_)+")");
		Vector<double> param(1,phi_);

		w.add_header()->title(s,'<');
		w<<param;
		w.add_header()->add(system_info_.get());
	} else { w<<phi_<<" "; }
}

Matrix<double> TrianglePhi::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 1;
	tmp(1,0) = 0;
	tmp(0,1) = 0.5;
	tmp(1,1) = sqrt(3.0)/2.0;
	return tmp;
}
/*}*/

/*{method needed for checking*/
void TrianglePhi::display_results(){
	compute_H();
	draw_lattice();

	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title(RST::math("\\phi="+my::tostring(phi_)));
		std::string run_cmd("./mc -s:wf triangle-free");
		run_cmd += " -u:N " + my::tostring(N_);
		run_cmd += " -u:m " + my::tostring(m_);
		run_cmd += " -u:n " + my::tostring(n_);
		run_cmd += " -i:bc "+ my::tostring(bc_);
		run_cmd += " -d:phi "+ my::tostring(phi_);
		run_cmd += " -d -u:tmax 10";

		rst_file_->title(title,'-');
		rst_file_->change_text_onclick("run command",run_cmd);

		rst_file_->figure(dir_+filename_+".png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(dir_+filename_+".pdf")+RST::scale("200"));
	}
}

void TrianglePhi::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="triangle-phi";
	display_results();
}
/*}*/
