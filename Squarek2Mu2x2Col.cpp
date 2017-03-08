#include "Squarek2Mu2x2Col.hpp"

Squarek2Mu2x2Col::Squarek2Mu2x2Col(System const& s, double const& mu, Vector<double> const& t):
	System(s),
	Square<double>(set_ab(),4,"square-k2mu-2x2-col"),
	mu_(mu),
	t_(t)
{
	if(N_/m_==2 && t_.size() == 4){
		if(status_==3){ init_lattice(); }
		if(status_==2){
			init_fermionic();
			same_wf_ = false;

			system_info_.text("Squarek2Mu2x2Col :");
			system_info_.item("Each color has a different Hamiltonian.");

			filename_ += "-mu"+std::string(mu_>=0?"+":"")+my::tostring(mu_);
			filename_ += "-t";
			for(unsigned int i(0);i<t_.size();i++){
				filename_ += ((t_(i)>=0)?"+":"")+my::tostring(t_(i));
			}
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : N/m != 2 or t_ wrong size : t_.size()="<<t_.size()<<std::endl; }
}

/*{method needed for running*/
void Squarek2Mu2x2Col::compute_H(unsigned int const& c){
	H_.set(n_,n_,0);

	unsigned int s0(0);
	unsigned int s1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		if(obs_[0](i,3)==1){
			if(obs_[0](i,5) < 2){ H_(s0,s1) = (obs_[0](i,4)?bc_*t_(1):t_(1)); }
			else { H_(s0,s1) = (obs_[0](i,4)?bc_*t_(3):t_(3)); }
		} else {
			switch(obs_[0](i,5)){
			case 0: { H_(s0,s1) =-(obs_[0](i,4)?bc_*t_(0):t_(0)); }break;
			case 1: { H_(s0,s1) = (obs_[0](i,4)?bc_*t_(2):t_(2)); }break;
			case 2: { H_(s0,s1) = (obs_[0](i,4)?bc_*t_(2):t_(2)); }break;
			case 3: { H_(s0,s1) = (obs_[0](i,4)?bc_*t_(0):t_(0)); }break;
			}
		}
		if((unsigned int)(obs_[0](i,5))%2==c%2){ H_(s0,s0) = mu_/2; }
	}
	H_ += H_.transpose();
}

void Squarek2Mu2x2Col::create(){
	for(unsigned int c(0);c<N_;c++){
		status_ = 2;
		compute_H(c);
		diagonalize(true);
		if(status_==1){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = H_(i,j);
				}
			}
		}
	}
}

void Squarek2Mu2x2Col::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("mu="+my::tostring(mu_)+" t=(");
		Vector<double> param(5);
		param(0) = mu_;
		for(unsigned int i(0);i<t_.size();i++){
			param(i+1) = t_(i);
			s += my::tostring(t_(i))+(i+1!=t_.size()?",":")");
		}

		w.add_to_header()->title(s,'<');
		w<<param;
		w.add_to_header()->add(system_info_.get());
	} else { w<<mu_<<" "<<t_; }
}

Matrix<double> Squarek2Mu2x2Col::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 2.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 0.0;
	tmp(1,1) = 2.0;
	return tmp;
}

unsigned int Squarek2Mu2x2Col::unit_cell_index(Vector<double> const& x) const {
	if(my::are_equal(x(0),0.0,eq_prec_,eq_prec_) && my::are_equal(x(1),0.0,eq_prec_,eq_prec_)){ return 0; }
	if(my::are_equal(x(0),0.5,eq_prec_,eq_prec_) && my::are_equal(x(1),0.0,eq_prec_,eq_prec_)){ return 1; }
	if(my::are_equal(x(0),0.5,eq_prec_,eq_prec_) && my::are_equal(x(1),0.5,eq_prec_,eq_prec_)){ return 2; }
	if(my::are_equal(x(0),0.0,eq_prec_,eq_prec_) && my::are_equal(x(1),0.5,eq_prec_,eq_prec_)){ return 3; }
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void Squarek2Mu2x2Col::display_results(){
	compute_H(0);
	draw_lattice(false,true,(spuc_==3?dir_nn_[2]*0.25+dir_nn_[3]*0.5:dir_nn_[2]*0.5+dir_nn_[3]*0.5));

	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title(RST::math("\\mu")+"="+my::tostring(mu_) + " ");
		std::string run_cmd("./mc -s:wf square-k2mu-2x2-col");
		run_cmd += " -u:N " + my::tostring(N_);
		run_cmd += " -u:m " + my::tostring(m_);
		run_cmd += " -u:n " + my::tostring(n_);
		run_cmd += " -i:bc "+ my::tostring(bc_);
		run_cmd += " -d:mu "+ my::tostring(mu_);
		run_cmd += " -d:t ";
		title += RST::math("t") + "=(";
		for(unsigned int i(0);i<t_.size();i++){
			run_cmd += my::tostring(t_(i))+(i+1!=t_.size()?",":"");
			title += my::tostring(t_(i))+(i+1!=t_.size()?",":")");
		}
		run_cmd += " -d -u:tmax 10";

		rst_file_->title(title,'-');
		rst_file_->change_text_onclick("run command",run_cmd);

		rst_file_->figure(dir_+filename_+".png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(dir_+filename_+".pdf")+RST::scale("200"));
	}
}

void Squarek2Mu2x2Col::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="square-k2mu-2x2-col";
	display_results();

	//plot_band_structure();
}
/*}*/
