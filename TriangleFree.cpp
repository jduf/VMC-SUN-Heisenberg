#include "TriangleFree.hpp"

TriangleFree::TriangleFree(System const& s, Vector<double> const& t, Vector<double> const& mu):
	System(s),
	Triangle<double>(set_ab(),3,"triangle-free"),
	t_(t),
	mu_(mu)
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		system_info_.text("TriangleFree :");
		system_info_.item("Each color has a the same Hamiltonian.");

		filename_ += "-t";
		for(unsigned int i(0);i<t_.size();i++){
			filename_ += ((t_(i)>=0)?"+":"")+my::tostring(t_(i));
		}
		filename_ += "-mu";
		for(unsigned int i(0);i<mu_.size();i++){
			filename_ += ((mu_(i)>=0)?"+":"")+my::tostring(mu_(i));
		}
	}
}

/*{method needed for running*/
void TriangleFree::compute_H(){
	H_.set(n_,n_,0);

	unsigned int s0(0);
	unsigned int s1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);

		switch(obs_[0](i,5)){
			case 0:
				{
					if(obs_[0](i,3)==2){ H_(s0,s1) = (obs_[0](i,4)?bc_:1)*t_(0); }
					else { H_(s0,s1) = (obs_[0](i,4)?bc_:1)*t_(1); }
					H_(s0,s0) = mu_(0);
				}break;
			case 1:
				{
					if(obs_[0](i,3)!=2){ H_(s0,s1) = (obs_[0](i,4)?bc_:1)*t_(0); }
					else { H_(s0,s1) = (obs_[0](i,4)?bc_:1)*t_(1); }
					H_(s0,s0) = mu_(1);
				}break;
			case 2:
				{
					H_(s0,s1) = (obs_[0](i,4)?bc_:1)*t_(0);
					H_(s0,s0) = mu_(2);
				}break;
		}
	}
	H_ += H_.transpose();
}

void TriangleFree::create(){
	compute_H();
	diagonalize(true);
	for(unsigned int c(0);c<N_;c++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_(c);j++){
				EVec_[c](i,j) = H_(i,j);
			}
		}
	}
}

void TriangleFree::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("t=(");
		Vector<double> param(t_.size()+mu_.size());

		for(unsigned int i(0);i<t_.size()-1;i++){
			param(i) = t_(i);
			s += my::tostring(t_(i))+",";
		}
		param(t_.size()-1) = t_.back();
		s += my::tostring(t_.back())+") "+RST::math("\\mu")+"=(";

		for(unsigned int i(0);i<mu_.size()-1;i++){
			param(i+t_.size()) = mu_(i);
			s += my::tostring(mu_(i))+",";
		}
		param.back() = mu_.back();
		s += my::tostring(mu_.back())+")";

		w.add_to_header()->title(s,'<');
		w<<param;
		w.add_to_header()->add(system_info_.get());
	} else { w<<t_<<" "<<mu_<<" "; }
}

Matrix<double> TriangleFree::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 1.5;
	tmp(1,0) =-sqrt(3.0)/2;
	tmp(0,1) = 1.5;
	tmp(1,1) = sqrt(3.0)/2;
	return tmp;
}

unsigned int TriangleFree::unit_cell_index(Vector<double> const& x) const {
	Vector<double> match(2,0);
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 0; }
	match(0) = 1.0/3.0;
	match(1) = 1.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 1; }
	match(0) = 2.0/3.0;
	match(1) = 2.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 2; }
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return 3;
}
/*}*/

/*{method needed for checking*/
void TriangleFree::display_results(){
	compute_H();
	draw_lattice(false,true,false,(dir_nn_[4]+dir_nn_[3])*1.5);

	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title("t=(");
		std::string run_cmd("./mc -s:wf triangle-free");
		run_cmd += " -u:N " + my::tostring(N_);
		run_cmd += " -u:m " + my::tostring(m_);
		run_cmd += " -u:n " + my::tostring(n_);
		run_cmd += " -i:bc "+ my::tostring(bc_);
		run_cmd += " -d:t ";
		for(unsigned int i(0);i<t_.size()-1;i++){
			title   += my::tostring(t_(i)) + ",";
			run_cmd += my::tostring(t_(i)) + ",";
		}
		title   += my::tostring(t_.back()) + "), "+RST::math("\\mu")+"=(";
		run_cmd += my::tostring(t_.back()) + " -d:mu ";
		for(unsigned int i(0);i<mu_.size()-1;i++){
			title   += my::tostring(mu_(i)) + ",";
			run_cmd += my::tostring(mu_(i)) + ",";
		}
		title   += my::tostring(mu_.back()) + ")";
		run_cmd += my::tostring(mu_.back()) + " -d -u:tmax 10";

		rst_file_->title(title,'-');
		rst_file_->change_text_onclick("run command",run_cmd);

		rst_file_->figure(dir_+filename_+".png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(dir_+filename_+".pdf")+RST::scale("200"));
	}
}

void TriangleFree::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="triangle-free";
	display_results();

	//compute_H();
	//plot_band_structure();
}
/*}*/
