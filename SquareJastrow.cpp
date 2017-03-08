#include "SquareJastrow.hpp"

SquareJastrow::SquareJastrow(System const& s, Matrix<double> const& nu):
	System(s),
	Square<double>(set_ab(),2,"square-jastrow")
{
	if(status_==3){ init_lattice(); }
	init_bosonic(z_,nu);

	system_info_.text("Staggered magnetic field, Becca's idea to mimic an on site chemical potential");
	system_info_.item("only works on a two sites per unit cell system (therefore only for SU(2))");

	status_=2;
}

/*{method needed for running*/
void SquareJastrow::create(){
	unsigned int l(z_/2);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		nn_(obs_[0](i,0),obs_[0](i,3)) = obs_[0](i,1);
		nn_(obs_[0](i,1),obs_[0](i,3)+l) = obs_[0](i,0);
		if(!(i%l)){ sl_(obs_[0](i,0)) = obs_[0](i,5); }
	}
	std::cout<<"the matrix should contain all neighbour for each site"<<std::endl;
	std::cout<<nn_<<std::endl;
	std::cout<<"the vector should contain the sublattice that the site belongs to"<<std::endl;
	std::cout<<sl_<<std::endl;

	if(N_==2){
		omega_(1,1) = -1.0;
		cc_(0,0) = 0;//|up,up>
		cc_(0,1) = 1;//|up,down>
		cc_(1,0) = 1;//|down,up>
		cc_(1,1) = 1;//|down,down>
	}
	/*!\warning omega might need to be complex*/
	//if(N_==3){
	//omega_(1,1) = std::polar(1.0,2.0*M_PI/3.0);
	//omega_(2,2) = std::polar(1.0,2.0*M_PI/3.0);
	//omega_(1,2) = std::polar(1.0,4.0*M_PI/3.0);
	//omega_(2,1) = std::polar(1.0,4.0*M_PI/3.0);
	//cc_(0,0) = 0;
	//cc_(0,1) = 1;
	//cc_(0,2) = 2;
	//cc_(1,0) = 1;
	//cc_(1,1) = 3;
	//cc_(1,2) = 4;
	//cc_(2,0) = 2;
	//cc_(2,1) = 4;
	//cc_(2,2) = 4;
	//}
	status_ = 1;
}

void SquareJastrow::save_param(IOFiles& w) const {
	w.write("nn (nearst neighbours)",nn_);
	w.write("cc (to match nu and x)",cc_);
	w.write("sl (sublattice)",sl_);
	w.write("omega (omega)",omega_);

	w.add_to_header()->add(system_info_.get());
}

Matrix<double> SquareJastrow::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 2;
	tmp(1,0) = 0;
	tmp(0,1) = 0;
	tmp(1,1) = 1;
	return tmp;
}

unsigned int SquareJastrow::unit_cell_index(Vector<double> const& x) const {
	return my::are_equal(x(0),0.5);
}
/*}*/

/*{method needed for checking*/
void SquareJastrow::display_results(){
	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title("Jastrow");
		std::string run_cmd("./mc -s:wf square-jastow");
		run_cmd += " -u:N " + my::tostring(N_);
		run_cmd += " -u:m " + my::tostring(m_);
		run_cmd += " -u:n " + my::tostring(n_);
		run_cmd += " -i:bc "+ my::tostring(bc_);
		run_cmd += " -d -u:tmax 10";

		rst_file_->title(title,'-');
		rst_file_->change_text_onclick("run command",run_cmd);

		rst_file_->figure(dir_+filename_+".png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(dir_+filename_+".pdf")+RST::scale("200"));
	}
}

void SquareJastrow::check(){
	std::cout<<"void SquareJastrow::check() : nothing to do"<<std::endl;
}
/*}*/
