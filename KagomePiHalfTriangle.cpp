#include "KagomePiHalfTriangle.hpp"

KagomePiHalfTriangle::KagomePiHalfTriangle(System const& s, Vector<double> const& t, Vector<double> const& phi):
	System(s),
	Kagome<std::complex<double> >(set_ab(),6,"kagome-pihalftriangle"),
	t_(t),
	phi_(phi)
{
	if(t_.size()==6){
		if(status_==3){ init_lattice(); }
		if(status_==2){
			init_fermionic();

			system_info_.text("KagomePiHalfTriangle :");
			system_info_.item("Each color has the same Hamiltonian.");
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : t must contain 6 values (currently contains "<<t_.size()<<")"<<std::endl; }
}

/*{method needed for running*/
void KagomePiHalfTriangle::compute_H(){
	H_.set(n_,n_,0);

	unsigned int t(0);
	unsigned int p(0);
	int s(0);
	for(unsigned int i(0);i<obs_[0].nlinks(); i++){
		switch(obs_[0](i,5)){
			case 0: { t = (obs_[0](i,3)==1?1:0); s = (obs_[0](i,4)?bc_:1); p = (obs_[0](i,3)==1?1:0); } break;
			case 1: { t = (obs_[0](i,3)==2?3:2); s = (obs_[0](i,4)?bc_:1); p = (obs_[0](i,3)==1?2:0); } break;
			case 2: { t = (obs_[0](i,3)==1?5:4); s = (obs_[0](i,4)?bc_:1); p = 0; } break;
			case 3: { t = (obs_[0](i,3)==1?1:0); s = (obs_[0](i,4)?bc_:1); p = (obs_[0](i,3)==1?2:0); } break;
			case 4: { t = (obs_[0](i,3)==2?3:2); s = (obs_[0](i,4)?bc_:1); p = (obs_[0](i,3)==1?2:0); } break;
			case 5: { t = (obs_[0](i,3)==1?5:4); s = (obs_[0](i,4)?bc_:1); p = (obs_[0](i,3)==2?3:0); } break;
		}
		H_(obs_[0](i,0),obs_[0](i,1)) = std::polar(s*t_(t),phi_(p)*M_PI);
	}
	H_ += H_.conjugate_transpose();
}

void KagomePiHalfTriangle::create(){
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

void KagomePiHalfTriangle::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("t=(");
		Vector<double> param(t_.size()+phi_.size());

		for(unsigned int i(0);i<t_.size()-1;i++){
			param(i) = t_(i);
			s += my::tostring(t_(i))+",";
		}
		param(t_.size()-1) = t_.back();
		s += my::tostring(t_.back())+") "+RST::math("\\phi")+"=(";

		for(unsigned int i(0);i<phi_.size()-1;i++){
			param(i+t_.size()) = phi_(i);
			s   += my::tostring(phi_(i)) + ",";
		}
		param.back() = phi_.back();
		s += my::tostring(phi_.back())+")";

		w.add_header()->title(s,'<');
		w<<param;
		w.add_header()->add(system_info_.get());
	} else { w<<t_<<" "; }
}

Matrix<double> KagomePiHalfTriangle::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 4.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 1.0;
	tmp(1,1) = sqrt(3.0);
	return tmp;
}

unsigned int KagomePiHalfTriangle::unit_cell_index(Vector<double> const& x) const {
	if(my::are_equal(x(1),0.0,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.00,eq_prec_,eq_prec_)){ return 0; }
		if(my::are_equal(x(0),0.25,eq_prec_,eq_prec_)){ return 2; }
		if(my::are_equal(x(0),0.50,eq_prec_,eq_prec_)){ return 3; }
		if(my::are_equal(x(0),0.75,eq_prec_,eq_prec_)){ return 5; }
	} else {
		if(my::are_equal(x(0),0.0,eq_prec_,eq_prec_)){ return 1; }
		if(my::are_equal(x(0),0.5,eq_prec_,eq_prec_)){ return 4; }
	}
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void KagomePiHalfTriangle::display_results(){
	compute_H();
	draw_lattice(false,true,(dir_nn_[2]+dir_nn_[3])*0.5);

	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title("t=(");
		std::string run_cmd("./mc -s:wf kagome-pihalftriangle");
		run_cmd += " -u:N " + my::tostring(N_);
		run_cmd += " -u:m " + my::tostring(m_);
		run_cmd += " -u:n " + my::tostring(n_);
		run_cmd += " -i:bc "+ my::tostring(bc_);
		run_cmd += " -d:t ";
		for(unsigned int i(0);i<t_.size()-1;i++){
			title   += my::tostring(t_(i)) + ",";
			run_cmd += my::tostring(t_(i)) + ",";
		}
		title   += my::tostring(t_.back()) + ") "+RST::math("\\phi")+"=(";
		run_cmd += my::tostring(t_.back()) + " -d:phi ";
		for(unsigned int i(0);i<phi_.size()-1;i++){
			title   += my::tostring(phi_(i)) + ",";
			run_cmd += my::tostring(phi_(i)) + ",";
		}
		title   += my::tostring(phi_.back()) + ")";
		run_cmd += my::tostring(phi_.back()) + " -d -u:tmax 10";

		rst_file_->title(title,'-');
		rst_file_->change_text_onclick("run command",run_cmd);

		rst_file_->figure(dir_+filename_+".png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(dir_+filename_+".pdf")+RST::scale("200"));
	}
}

void KagomePiHalfTriangle::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="kagome-pihalftriangle";
	display_results();
}
/*}*/
