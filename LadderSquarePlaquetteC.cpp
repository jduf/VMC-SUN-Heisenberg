#include "LadderSquarePlaquetteC.hpp"

LadderSquarePlaquetteC::LadderSquarePlaquetteC(System const& s, Vector<double> const& t):
	System(s),
	Ladder<double>(4,"ladder-squareplaquetteC"),
	t_(t)
{
	if(status_==2 && t_.size()==3){
		init_fermionic();

		system_info_.text("LadderSquarePlaquetteC :");
		system_info_.text(" Each color has the same Hamiltonian.");
		system_info_.text(" Square plaquette in a 4 site unit cell");
		system_info_.text(" pi flux between the plaquettes");
		system_info_.text(" pi flux inside the plaquettes");

		filename_ += "-t";
		for(unsigned int i(0);i<t_.size();i++){
			filename_ += ((t_(i)>=0)?"+":"")+my::tostring(t_(i));
		}
	}
}

/*{method needed for running*/
void LadderSquarePlaquetteC::compute_H(){
	H_.set(n_,n_,0);

	double t(0.0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		switch(obs_[0](i,5)){
			case 0:
				{
					if(obs_[0](i,3)){ t = t_(1); }
					else            { t = t_(0); }
				}break;
			case 1:
				{ t = -t_(0); }break;
			case 2:
				{ 
					if(obs_[0](i,3)){ t = t_(2); }
					else            { t = t_(0); }
				}break;
			case 3:
				{ t = -t_(0); }break;
		}

		H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*t:t);
	}
	H_ += H_.transpose();
}

void LadderSquarePlaquetteC::create(){
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

void LadderSquarePlaquetteC::save_param(IOFiles& w) const {
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
/*}*/

/*{method needed for checking*/
void LadderSquarePlaquetteC::display_results(){
	compute_H();
	draw_lattice(true,true,true);

	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title(RST::math("\\theta=")+my::tostring(acos(J_(0))) + " : t=(");
		std::string run_cmd("./mc -s:wf ladder-squareplaquetteC");
		run_cmd += " -u:N " + my::tostring(N_);
		run_cmd += " -u:m " + my::tostring(m_);
		run_cmd += " -u:n " + my::tostring(n_);
		run_cmd += " -i:bc "+ my::tostring(bc_);
		run_cmd += " -d:theta " + my::tostring(acos(J_(0)));
		run_cmd += " -d:t ";
		for(unsigned int i(0);i<t_.size()-1;i++){
			title   += my::tostring(t_(i)) + ",";
			run_cmd += my::tostring(t_(i)) + ",";
		}
		title   += my::tostring(t_.back()) + ")";
		run_cmd += my::tostring(t_.back()) + " -d -u:tmax 10";
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

void LadderSquarePlaquetteC::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="ladder-squareplaquetteC";
	display_results();

	//compute_H();
	//plot_band_structure();
}
/*}*/
