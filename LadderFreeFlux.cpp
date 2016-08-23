#include"LadderFreeFlux.hpp"

LadderFreeFlux::LadderFreeFlux(System const& s, Vector<double> const& t, Vector<double> const& flux):
	System(s),
	Ladder<std::complex<double> >(set_spuc(t,flux),"ladder-freeflux"),
	t_(t),
	flux_(flux)
{
	if(status_==2){
		init_fermionic();

		system_info_.text("LadderFreeFlux :");
		system_info_.item("Each color has the same Hamiltonian.");
		system_info_.item("Free hopping term");
		system_info_.item("Free flux per plaquette ");
		filename_ += "-t";
		for(unsigned int i(0);i<t_.size();i++){
			filename_ += ((t_(i)>0)?"+":"")+my::tostring(t_(i));
		}
		filename_ += "-phi";
		for(unsigned int i(0);i<flux_.size()-1;i++){
			filename_ += ((flux_(i)>0)?"+":"")+my::tostring(flux_(i));
		}
	}
}

/*{method needed for running*/
void LadderFreeFlux::compute_H(){
	H_.set(n_,n_,0);
	Matrix<int> nb;
	std::cerr<<__PRETTY_FUNCTION__<<" : undefined"<<std::endl;
	//for(unsigned int i(0);i<n_;i++){
		//nb = get_neighbourg(i);
		//if(i%spuc_){
			//if(i%2){
				//H_(i,nb(0,0)) = my::chop(std::polar(nb(0,1)*t_(k++),flux_(f++)*M_PI));
			//} else {
				//H_(i,nb(0,0)) = nb(0,1)*t_(k++);
				//H_(i,nb(1,0)) = nb(1,1)*t_(k++);
			//}
		//} else {
			///*!decoupled chains in the limit J⊥(0)=J_(1)->0, therefore, in
			// * that case the variational parameter is t⊥ (t‖ is set to 1),
			// * otherwise the inverse is done*/
			//if(J_(0)>J_(1)){
				//H_(i,nb(0,0)) = nb(0,1);
				//H_(i,nb(1,0)) = nb(1,1)*t_(k++);
			//} else {
				//H_(i,nb(0,0)) = nb(0,1)*t_(k++);
				//H_(i,nb(1,0)) = nb(1,1);
			//}
		//}
		//k = k%t_.size();
		//f = f%flux_.size();
	//}
	H_ += H_.conjugate_transpose();
}

void LadderFreeFlux::create(){
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

void LadderFreeFlux::save_param(IOFiles& w) const {
	if(w.is_binary()){
		Vector<double> param(t_.size()+flux_.size());
		std::string s("t=(");
		for(unsigned int i(0);i<t_.size()-1;i++){
			s += my::tostring(t_(i))+",";
			param(i) = t_(i);
		}
		param(t_.size()-1) = t_.back();
		s += my::tostring(t_.back())+") : " + RST::math("\\phi/\\pi")+"=(";
		for(unsigned int i(0);i<flux_.size()-1;i++){
			s += my::tostring(flux_(i))+",";
			param(i+t_.size()) = flux_(i);
		}
		param.back() = flux_.back();
		s += my::tostring(flux_.back())+")";
		w.add_header()->title(s,'<');
		w<<param;
		w.add_header()->add(system_info_.get());
	} else { w<<t_<<" "<<flux_<<" "; }
}

unsigned int LadderFreeFlux::set_spuc(Vector<double> const& t, Vector<double> const& flux){
	if((t.size()+1)%flux.size()){
		std::cerr<<__PRETTY_FUNCTION__<<" : incoherent t and flux sizes : "<<t.size()<<" "<<flux.size()<<std::endl;
		return 1;
	} else {
		switch(t.size()){
			case 2: { return 2; } break;
			case 5: { return 4; } break;
			case 8: { return 6; } break;
			case 11:{ return 8; } break;
			default:{
						std::cerr<<__PRETTY_FUNCTION__<<" : invalid t size : "<<t.size()<<std::endl;
						return 1;
					}
		}
	}
}
/*}*/

/*{method needed for checking*/
void LadderFreeFlux::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="ladder-freeflux";
	display_results();

	//plot_band_structure();
}

void LadderFreeFlux::display_results(){
	compute_H();
	draw_lattice(true,true,true);

	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title(RST::math("\\theta=")+my::tostring(acos(J_(0))) + " : t=(");
		std::string run_cmd("./mc -s:wf ladder-free");
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
		title   += my::tostring(t_.back()) + "), "+RST::math("\\mu")+"=(";
		run_cmd += my::tostring(t_.back()) + " -d:phi ";
		for(unsigned int i(0);i<flux_.size()-1;i++){
			title   += my::tostring(flux_(i)) + ",";
			run_cmd += my::tostring(flux_(i)) + ",";
		}
		title   += my::tostring(flux_.back()) + ")";
		run_cmd += my::tostring(flux_.back()) + " -d -u:tmax 10";
		if(dir_ == "P/" || dir_ == "O/" || dir_ == "A/"){
			rst_file_->title("|theta"+my::tostring(acos(J_(0)))+"|_",'-');
			rst_file_->replace("theta"+my::tostring(acos(J_(0))),title);
		} else { rst_file_->title(title,'-'); }

		rst_file_->change_text_onclick("run command",run_cmd);
		rst_file_->figure(dir_+filename_+".png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(dir_+filename_+".pdf")+RST::scale("200"));

		Vector<unsigned int> o(6,0);
		unsigned int o_index_tmp(2);
		for(unsigned int i(1);i<obs_.size();i++){
			switch(obs_[i].get_type()){
				case 1:{ o(0)=i; }break;//bond energy
				case 2:{ o(o_index_tmp++)=i; }break;//long range correlation
				case 3:{ o(1)=i; }break;//color occupation
			}
		}
		if(o(2) && o(3)){ rst_file_->figure(relative_path+filename_+"-lr.png","long range correlations",RST::target(relative_path+filename_+"-lr.gp")+RST::scale("200")); }
		if(o(2) && o(3) && o(4) && o(5)){ rst_file_->figure(relative_path+filename_+"-as.png","(anti)symmetric correlations",RST::target(relative_path+filename_+"-as.gp")+RST::scale("200")); }
	}
}
/*}*/

/*{method needed for analysing*/
std::string LadderFreeFlux::extract_level_6(){
	(*data_write_)<<N_<<" "<<m_<<" "<<n_<<" "<<bc_<<" "<<asin(J_(1))<<" "<<obs_[0][0]<<IOFiles::endl;

	display_results();

	save_param(*jd_write_);
	save(*jd_write_);

	return filename_;
}
/*}*/
