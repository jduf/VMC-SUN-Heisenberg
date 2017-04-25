#include "TriangleAlternatingPlaquette.hpp"

TriangleAlternatingPlaquette::TriangleAlternatingPlaquette(System const& s, double const& t):
	System(s),
	Triangle<double>(set_ab(),6,"triangle-alternatingplaquette"),
	t_(t)
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		system_info_.text("TriangleAlternatingPlaquette :");
		system_info_.item("Each color has the same Hamiltonian.");
		system_info_.item("Triangular plaquette with different hopping term than the rest of the lattice");

		filename_ += "-t"+std::string(t_>=0?"+":"")+my::tostring(t_);
	}
}

/*{method needed for running*/
void TriangleAlternatingPlaquette::compute_H(){
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
					if(obs_[0](i,3)!=2){ H_(s0,s1) = (obs_[0](i,4)?bc_*t:t); }
					else { H_(s0,s1) = (obs_[0](i,4)?bc_*t_:t_); }
				}break;
			case 1:
				{
					if(obs_[0](i,3)!=2){ H_(s0,s1) = (obs_[0](i,4)?bc_*t_:t_); }
					else { H_(s0,s1) = (obs_[0](i,4)?bc_*t:t); }
				}break;
			case 2:
				{
					if(obs_[0](i,3)!=0){ H_(s0,s1) = (obs_[0](i,4)?bc_*t:t); }
					else { H_(s0,s1) = (obs_[0](i,4)?bc_*t_:t_); }
				}break;
			case 4:
				{
					if(obs_[0](i,3)!=0){ H_(s0,s1) = (obs_[0](i,4)?bc_*t_:t_); }
					else { H_(s0,s1) = (obs_[0](i,4)?bc_*t:t); }
				}break;
			case 3:
			case 5:
				{ H_(s0,s1) = (obs_[0](i,4)?bc_*t_:t_); }break;
		}
	}
	H_ += H_.transpose();
}

void TriangleAlternatingPlaquette::create(){
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

void TriangleAlternatingPlaquette::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("t="+my::tostring(t_));
		Vector<double> param(1,t_);

		w.add_to_header()->title(s,'<');
		w<<param;
		w.add_to_header()->add(system_info_.get());
	} else { w<<t_<<" "; }
}

Matrix<double> TriangleAlternatingPlaquette::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 3.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 0.0;
	tmp(1,1) = sqrt(3.0);
	return tmp;
}

unsigned int TriangleAlternatingPlaquette::unit_cell_index(Vector<double> const& x) const {
	if(my::are_equal(x(1),0.0,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.0,    eq_prec_,eq_prec_)){ return 0; }
		if(my::are_equal(x(0),1.0/3.0,eq_prec_,eq_prec_)){ return 1; }
		if(my::are_equal(x(0),2.0/3.0,eq_prec_,eq_prec_)){ return 2; }
	}
	if(my::are_equal(x(1),0.5,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),1.0/6.0,eq_prec_,eq_prec_)){ return 3; }
		if(my::are_equal(x(0),3.0/6.0,eq_prec_,eq_prec_)){ return 4; }
		if(my::are_equal(x(0),5.0/6.0,eq_prec_,eq_prec_)){ return 5; }
	}
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void TriangleAlternatingPlaquette::display_results(){
	compute_H();
	draw_lattice(true,true,false,(dir_nn_[3]+dir_nn_[4])*0.25);

	if(rst_file_){
		std::string title("t="+my::tostring(t_));
		rst_file_->title(title,'-');
		rst_file_->change_text_onclick("run command",get_mc_run_command());

		rst_file_->figure(dir_+filename_+".png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(dir_+filename_+".pdf")+RST::scale("200"));
	}
}

void TriangleAlternatingPlaquette::check(){
	//info_ = "";
	//path_ = "";
	//dir_  = "./";
	//filename_ ="triangle-alternatingplaquette";
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

std::string TriangleAlternatingPlaquette::get_mc_run_command() const {
	std::string run_cmd("./mc -s:wf triangle-alternatingplaquette");
	run_cmd += " -u:N " + my::tostring(N_);
	run_cmd += " -u:m " + my::tostring(m_);
	run_cmd += " -u:n " + my::tostring(n_);
	run_cmd += " -i:bc "+ my::tostring(bc_);
	run_cmd += " -d:t " + my::tostring(t_);
	run_cmd += " -d -u:tmax 10";

	return run_cmd;
}
/*}*/
