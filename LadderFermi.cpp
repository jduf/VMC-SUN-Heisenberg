#include"LadderFermi.hpp"

LadderFermi::LadderFermi(System const& s):
	System(s),
	Ladder(1,"ladder-fermi")
{
	if(status_==2){
		init_fermionic();

		system_info_.text("LadderFermi :");
		system_info_.item("Each color has the same Hamiltonian.");
		system_info_.item("Uniform real hopping term.");
	}
}

/*{method needed for running*/
void LadderFermi::compute_H(){
	H_.set(n_,n_,0);
	for(unsigned int i(0);i<obs_[0].nlinks(); i++){
		H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_:1);
	}
	H_ += H_.transpose();
}

void LadderFermi::create(){
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
/*}*/

/*{method needed for checking*/
void LadderFermi::display_results(){
	compute_H();
	draw_lattice(true,true,true);

	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title(RST::math("\\theta=")+my::tostring(acos(J_(0))) + " : Fermi");
		std::string run_cmd("./mc -s:wf ladder-fermi");
		run_cmd += " -u:N " + my::tostring(N_);
		run_cmd += " -u:m " + my::tostring(m_);
		run_cmd += " -u:n " + my::tostring(n_);
		run_cmd += " -i:bc "+ my::tostring(bc_);
		run_cmd += " -d:theta " + my::tostring(acos(J_(0)));
		run_cmd += " -d -u:tmax 10";
		if(dir_ == "P/" || dir_ == "O/" || dir_ == "A/"){
			rst_file_->title("|theta"+my::tostring(acos(J_(0)))+"|_",'-');
			rst_file_->replace("theta"+my::tostring(acos(J_(0))),title);
		} else { rst_file_->title(title,'-'); }

		rst_file_->change_text_onclick("run command",run_cmd);
		rst_file_->figure(dir_+filename_+".png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(dir_+filename_+".pdf")+RST::scale("200"));
		rst_file_->figure(relative_path+filename_+"-lr.png","long range correlations",RST::target(relative_path+filename_+"-lr.gp")+RST::scale("200")); 
		rst_file_->figure(relative_path+filename_+"-as.png","(anti)symmetric correlations",RST::target(relative_path+filename_+"-as.gp")+RST::scale("200")); 
	}
}

void LadderFermi::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="ladder-fermi";
	display_results();

	//plot_band_structure();
}
/*}*/
