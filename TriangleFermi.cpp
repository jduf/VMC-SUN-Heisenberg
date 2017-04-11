#include "TriangleFermi.hpp"

TriangleFermi::TriangleFermi(System const& s):
	System(s),
	Triangle<double>(set_ab(),1,"triangle-fermi")
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		system_info_.text("TriangleFermi :");
		system_info_.item("Each color has the same Hamiltonian.");
	}
}

/*{method needed for running*/
void TriangleFermi::compute_H(){
	H_.set(n_,n_,0);

	double t(-1.0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*t:t);
	}
	H_ += H_.transpose();
}

void TriangleFermi::create(){
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

Matrix<double> TriangleFermi::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 1.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 0.5;
	tmp(1,1) = sqrt(3.0)/2.0;
	return tmp;
}
/*}*/

/*{method needed for checking*/
void TriangleFermi::display_results(){
	compute_H();
	draw_lattice(false,true,false,(dir_nn_[4]+dir_nn_[3])*1.5);

	if(rst_file_){
		std::string title("Fermi");
		std::string run_cmd("./mc -s:wf triangle-fermi");
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

void TriangleFermi::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="triangle-fermi";
	display_results();

	//plot_band_structure();
	//Vector<unsigned int> c(3,0);
	//for(unsigned int i(0);i<obs_[0].nlinks();i++){
		//c(obs_[0](i,2))++;
		//if(obs_[0](i,2)==2){ std::cout<<obs_[0](i,0)<<" "<<obs_[0](i,1)<<std::endl; }
	//}
	//create_obs(2);
	//Vector<unsigned int> c(n_,0);
	//for(auto const& o:obs_){
		//if(o.get_type()==2){
			//for(unsigned int i(0);i<o.nlinks();i++){
				//c(o(i,2))++;
			//}
		//}
	//}
	//std::cout<<c<<std::endl;
}
/*}*/
