#include "SquareT4x4.hpp"

SquareT4x4::SquareT4x4(System const& s, Vector<double> const& t):
	System(s),
	Square<double>(set_ab(),16,"square-T4x4"),
	t_(t)
{
	if(t_.size()==32){
		if(status_==3){ init_lattice(); }
		if(status_==2){
			init_fermionic();

			system_info_.text("SquareT4x4 :");
			system_info_.item("Each color has the same Hamiltonian.");
			system_info_.item("16 sites in a 4x4 unit cell");

			filename_ += "-t";
			for(unsigned int i(0);i<t_.size();i++){
				filename_ += ((t_(i)>=0)?"+":"")+my::tostring(t_(i));
			}
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : t must contain 32 values (currently contains "<<t_.size()<<")"<<std::endl; }
}

/*{method needed for running*/
void SquareT4x4::compute_H(){
	H_.set(n_,n_,0);

	double t(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		switch(obs_[0](i,5)){
			case 0: { t = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)?1:0); }break;
			case 1: { t = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)?3:2); }break;
			case 2: { t = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)?5:4); }break;
			case 3: { t = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)?7:6); }break;
			case 4: { t = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)?9:8); }break;
			case 5: { t = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)?11:10); }break;
			case 6: { t = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)?13:12); }break;
			case 7: { t = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)?15:14); }break;
			case 8: { t = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)?17:16); }break;
			case 9: { t = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)?19:18); }break;
			case 10:{ t = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)?21:20); }break;
			case 11:{ t = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)?23:22); }break;
			case 12:{ t = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)?25:24); }break;
			case 13:{ t = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)?27:26); }break;
			case 14:{ t = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)?29:28); }break;
			case 15:{ t = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)?31:30); }break;
		}
		H_(obs_[0](i,0),obs_[0](i,1)) = t;
	}
	H_ += H_.transpose();
}

void SquareT4x4::create(){
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

void SquareT4x4::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("t=(");
		Vector<double> param(t_.size());

		for(unsigned int i(0);i<t_.size()-1;i++){
			param(i) = t_(i);
			s += my::tostring(t_(i))+",";
		}
		param(t_.size()-1) = t_.back();
		s += my::tostring(t_.back())+")";

		w.add_to_header()->title(s,'<');
		w<<param;
		w.add_to_header()->add(system_info_.get());
	} else { w<<t_<<" "; }
}

Matrix<double> SquareT4x4::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 4.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 0.0;
	tmp(1,1) = 4.0;
	return tmp;
}

unsigned int SquareT4x4::unit_cell_index(Vector<double> const& x) const {
	if(my::are_equal(x(1),0.00,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.00,eq_prec_,eq_prec_)){ return 0; }
		if(my::are_equal(x(0),0.25,eq_prec_,eq_prec_)){ return 1; }
		if(my::are_equal(x(0),0.50,eq_prec_,eq_prec_)){ return 2; }
		if(my::are_equal(x(0),0.75,eq_prec_,eq_prec_)){ return 3; }
	}
	if(my::are_equal(x(1),0.25,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.00,eq_prec_,eq_prec_)){ return 4; }
		if(my::are_equal(x(0),0.25,eq_prec_,eq_prec_)){ return 5; }
		if(my::are_equal(x(0),0.50,eq_prec_,eq_prec_)){ return 6; }
		if(my::are_equal(x(0),0.75,eq_prec_,eq_prec_)){ return 7; }
	}
	if(my::are_equal(x(1),0.50,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.00,eq_prec_,eq_prec_)){ return 8; }
		if(my::are_equal(x(0),0.25,eq_prec_,eq_prec_)){ return 9; }
		if(my::are_equal(x(0),0.50,eq_prec_,eq_prec_)){ return 10;}
		if(my::are_equal(x(0),0.75,eq_prec_,eq_prec_)){ return 11;}
	}
	if(my::are_equal(x(1),0.75,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.00,eq_prec_,eq_prec_)){ return 12;}
		if(my::are_equal(x(0),0.25,eq_prec_,eq_prec_)){ return 13;}
		if(my::are_equal(x(0),0.50,eq_prec_,eq_prec_)){ return 14;}
		if(my::are_equal(x(0),0.75,eq_prec_,eq_prec_)){ return 15;}
	}
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void SquareT4x4::display_results(){
	compute_H();
	draw_lattice(true,true,false,(dir_nn_[2]+dir_nn_[3])*0.5);

	if(rst_file_){
		std::string title("t=(");
		std::string run_cmd("./mc -s:wf square-T4x4");
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

void SquareT4x4::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="square-T4x4";
	display_results();

	//plot_band_structure();
}
/*}*/
