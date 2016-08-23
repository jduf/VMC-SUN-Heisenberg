#include "TrianglePlaquette.hpp"

TrianglePlaquette::TrianglePlaquette(System const& s, double const& t):
	System(s),
	Triangle<double>(set_ab(),3,"triangle-plaquette"),
	t_(t)
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		system_info_.text("TrianglePlaquette :");
		system_info_.item("Each color has the same Hamiltonian.");
		system_info_.item("Triangular plaquette with different hopping term than the rest of the lattice");

		filename_ += "-t"+std::string(t_>=0?"+":"")+my::tostring(t_);
	}
}

/*{method needed for running*/
void TrianglePlaquette::compute_H(){
	H_.set(n_,n_,0);

	double t(-1.0);
	unsigned int s0(0);
	unsigned int s1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		switch(obs_[0](i,5)){
			case 0:
				{
					if(obs_[0](i,3)==2){ H_(s0,s1) = (obs_[0](i,4)?bc_*t_:t_); }
					else { H_(s0,s1) = (obs_[0](i,4)?bc_*t:t); }
				}break;
			case 1:
				{
					if(obs_[0](i,3)!=2){ H_(s0,s1) = (obs_[0](i,4)?bc_*t_:t_); }
					else { H_(s0,s1) = (obs_[0](i,4)?bc_*t:t); }
				}break;
			case 2:
				{ H_(s0,s1) = (obs_[0](i,4)?bc_*t_:t_); }break;
		}
	}
	H_ += H_.transpose();
}

void TrianglePlaquette::create(){
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

void TrianglePlaquette::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("t="+my::tostring(t_));
		Vector<double> param(1,t_);

		w.add_header()->title(s,'<');
		w<<param;
		w.add_header()->add(system_info_.get());
	} else { w<<t_<<" "; }
}

Matrix<double> TrianglePlaquette::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 1.5;
	tmp(1,0) =-sqrt(3.0)/2;
	tmp(0,1) = 1.5;
	tmp(1,1) = sqrt(3.0)/2;
	return tmp;
}

unsigned int TrianglePlaquette::unit_cell_index(Vector<double> const& x) const {
	Vector<double> match(2,0);
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 0; }
	match(0) = 1.0/3.0;
	match(1) = 1.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 1; }
	match(0) = 2.0/3.0;
	match(1) = 2.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 2; }
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void TrianglePlaquette::display_results(){
	compute_H();
	draw_lattice(true,true,dir_nn_[3]*0.5);

	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title("t="+my::tostring(t_));
		rst_file_->title(title,'-');
		rst_file_->change_text_onclick("run command",get_mc_run_command());

		rst_file_->figure(dir_+filename_+".png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(dir_+filename_+".pdf")+RST::scale("200"));
	}
}

void TrianglePlaquette::check(){
	//info_ = "";
	//path_ = "";
	//dir_  = "./";
	//filename_ ="triangle-plaquette";
	//display_results();

	//compute_H();
	//plot_band_structure();
	
	Data<double> b1;
	Data<double> b2;
	for(unsigned int i(0);i<obs_[1].nval();i++){
		if(obs_[1][i].get_x()<0){ b1.merge(obs_[1][i]); }
		else                    { b2.merge(obs_[1][i]); }
	}
	b1.complete_analysis(1e-5);
	b2.complete_analysis(1e-5);
	
	std::cerr<<N_<<" "<<m_<<" "<<n_<<" "<<bc_<<" "<<t_<<" "<<obs_[0][0]<<" "<<b1<<" "<<b2<<" "<<ref_<<std::endl;
}

std::string TrianglePlaquette::get_mc_run_command() const {
	std::string run_cmd("./mc -s:wf triangle-plaquette");
	run_cmd += " -u:N " + my::tostring(N_);
	run_cmd += " -u:m " + my::tostring(m_);
	run_cmd += " -u:n " + my::tostring(n_);
	run_cmd += " -i:bc "+ my::tostring(bc_);
	run_cmd += " -d:t " + my::tostring(t_);
	run_cmd += " -d -u:tmax 10";

	return run_cmd;
}
/*}*/
