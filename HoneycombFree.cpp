#include "HoneycombFree.hpp"

HoneycombFree::HoneycombFree(System const& s, Vector<double> t):
	System(s),
	Honeycomb<double>(set_ab(),6,"honeycomb-free"),
	t_(t)
{
	if(status_==2){
		init_fermionic();

		filename_ += "-t";
		for(unsigned int i(0);i<t_.size();i++){
			filename_ += ((t_(i)>=0)?"+":"")+my::tostring(t_(i));
		}
	}
}

/*{method needed for running*/
void HoneycombFree::compute_H(){
	H_.set(n_,n_,0);

	unsigned int s0(0);
	unsigned int s1(0);
	unsigned int ab0(0);
	unsigned int ab1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		ab0 = get_site_in_ab(s0);
		ab1 = get_site_in_ab(s1);
		if(ab0==0 && ab1==1){ H_(s0,s1) = obs_[0](i,4)*t_(0); }
		if(ab0==0 && ab1==3){ H_(s0,s1) = obs_[0](i,4)*t_(1); }
		if(ab0==0 && ab1==5){ H_(s0,s1) = obs_[0](i,4)*t_(2); }
		if(ab0==2 && ab1==1){ H_(s0,s1) = obs_[0](i,4)*t_(3); }
		if(ab0==2 && ab1==3){ H_(s0,s1) = obs_[0](i,4)*t_(4); }
		if(ab0==2 && ab1==5){ H_(s0,s1) = obs_[0](i,4)*t_(5); }
		if(ab0==4 && ab1==1){ H_(s0,s1) = obs_[0](i,4)*t_(6); }
		if(ab0==4 && ab1==3){ H_(s0,s1) = obs_[0](i,4)*t_(7); }
		if(ab0==4 && ab1==5){ H_(s0,s1) = obs_[0](i,4)*t_(8); }
	}
	H_ += H_.transpose();
}

void HoneycombFree::create(){
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

void HoneycombFree::save_param(IOFiles& w) const {
	std::string s("t=(");

	for(unsigned int i(0);i<t_.size()-1;i++){ 
		s += my::tostring(t_(i))+",";
	}
	s += my::tostring(t_.back())+")";

	w.add_header()->title(s,'<');
	w<<t_;
	GenericSystem<double>::save_param(w);
}

unsigned int HoneycombFree::match_pos_in_ab(Vector<double> const& x) const {
	Vector<double> match(2,0);
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 0; }
	match(0) = 1.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 1; }
	match(1) = 1.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 2; }
	match(0) = 0;
	match(1) = 2.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 3; }
	match(0) = 2.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 4; }
	match(1) = 1.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 5; }
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return 6;
}

Matrix<double> HoneycombFree::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 3.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 1.5;
	tmp(1,1) = 1.5*sqrt(3.0);
	return tmp;
}
/*}*/

/*{method needed for checking*/
void HoneycombFree::lattice(){
	compute_H();

	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(info_+path_+dir_,filename_+"-pstricks");
	ps.begin(-2,-20,40,20,filename_+"-pstricks");

	Matrix<double> polygon(4,2);
	polygon(0,0)=0;
	polygon(0,1)=0;
	polygon(1,0)=LxLy_(0,0);
	polygon(1,1)=LxLy_(1,0);
	polygon(2,0)=LxLy_(0,0)+LxLy_(0,1);
	polygon(2,1)=LxLy_(1,0)+LxLy_(1,1);
	polygon(3,0)=LxLy_(0,1);
	polygon(3,1)=LxLy_(1,1);
	ps.polygon(polygon,"linecolor=green");

	polygon(0,0)=0;
	polygon(0,1)=0;
	polygon(1,0)=ab_(0,0);
	polygon(1,1)=ab_(1,0);
	polygon(2,0)=ab_(0,0)+ab_(0,1);
	polygon(2,1)=ab_(1,0)+ab_(1,1);
	polygon(3,0)=ab_(0,1);
	polygon(3,1)=ab_(1,1);
	ps.polygon(polygon,"linecolor=black");

	double x_shift((LxLy_(0,0)+LxLy_(0,1)-ab_(0,0)-ab_(0,1))/2);
	double y_shift((-ab_(1,0)-ab_(1,1))/2);
	polygon(0,0)+=x_shift;
	polygon(0,1)+=y_shift;
	polygon(1,0)+=x_shift;
	polygon(1,1)+=y_shift;
	polygon(2,0)+=x_shift;
	polygon(2,1)+=y_shift;
	polygon(3,0)+=x_shift;
	polygon(3,1)+=y_shift;
	ps.polygon(polygon,"linecolor=black");

	double t;
	double corr;
	unsigned int s0;
	unsigned int s1;
	std::string str;
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		xy0 = get_pos_in_lattice(s0);

		s1 = obs_[0](i,1);
		xy1 = get_pos_in_lattice(s1);

		//if((my::in_polygon(polygon.row(),polygon.ptr(),polygon.ptr()+polygon.row(),xy0(0),xy0(1)) || my::in_polygon(polygon.row(),polygon.ptr(),polygon.ptr()+polygon.row(),xy1(0),xy1(1))) ){
			t = H_(s0,s1);
			if(std::abs(t)>1e-4){
				if((xy0-xy1).norm_squared()>1.0001){
					linestyle = "dashed";
					xy1 = xy0;
					xy1(0) += dir_nn_(obs_[0](i,3),0);
					xy1(1) += dir_nn_(obs_[0](i,3),1);
					xy1 = xy1.chop();
					ps.put(xy1(0)-0.20,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
				} else { 
					linestyle = "solid";  
					if(s0<s1){
						ps.put(xy0(0)-0.20,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}"); 
						ps.put(xy1(0)-0.20,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}"); 
					}
				}

				if(t>0){ color = "blue";}
				else { color = "red"; }
				linewidth = my::tostring(std::abs(t))+"mm";
				ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
			}
		//} else {
			if(obs_.size()>1){/*bound energy*/
				corr = obs_[1][obs_[0](i,2)].get_x();
				if(std::abs(corr)>1e-4){
					if(corr<0){ color = "red"; }
					else { color = "blue"; }
					linewidth = my::tostring(std::abs(corr))+"mm";

					ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
				}

				//if(i%3!=1){
				//if(i%3==0){
				//ps.put(xy0(0)+x_shift,xy0(1)-0.2,"\\tiny{"+my::tostring(s0)+"}");
				//}
				//str = my::tostring(corr);
				//ps.put((xy0(0)+xy1(0))/2.0+2*x_shift,xy0(1),"\\tiny{"+str.substr(0,8)+"}");
				//str = my::tostring(obs_[1][i].get_dx());
				////if(obs_[1][i].get_dx()<1e-4){
				////ps.put((xy0(0)+xy1(0))/2.0+2*x_shift,xy0(1)-0.2,"\\tiny{"+str.substr(0,4)+"e-"+str.substr(str.size()-2,2)+"}");
				////} else {
				////ps.put((xy0(0)+xy1(0))/2.0+2*x_shift,xy0(1)-0.2,"\\tiny{"+str.substr(0,8)+"}");
				////}
				//} else {
				//ps.put(xy1(0)+x_shift,xy1(1)+0.2,"\\tiny{"+my::tostring(s1)+"}");
				//str = my::tostring(corr);
				//ps.put((xy0(0)+xy1(0))/2.0+2*x_shift,(xy0(1)+xy1(1))/2.0,"\\tiny{"+str.substr(0,8)+"}");
				//str = my::tostring(obs_[1][i].get_dx());
				////if(obs_[1][i].get_dx()<1e-4){
				////ps.put((xy0(0)+xy1(0))/2.0+2*x_shift,(xy0(1)+xy1(1))/2.0-0.2,"\\tiny{"+str.substr(0,4)+"e-"+str.substr(str.size()-2,2)+"}");
				////} else {
				////ps.put((xy0(0)+xy1(0))/2.0+2*x_shift,(xy0(1)+xy1(1))/2.0-0.2,"\\tiny{"+str.substr(0,8)+"}");
				////}
				//}
			}
		//}
	}
	ps.end(true,true,true);
}

void HoneycombFree::display_results(){
	lattice();
	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title("t=(");
		std::string run_cmd("./mc -s:wf ladder-free");
		run_cmd += " -u:N " + my::tostring(N_);
		run_cmd += " -u:m " + my::tostring(m_);
		run_cmd += " -u:n " + my::tostring(n_);
		run_cmd += " -i:bc "+ my::tostring(bc_);
		run_cmd += " -d:t ";
		for(unsigned int i(0);i<t_.size()-1;i++){
			title   += my::tostring(t_(i)) + ","; 
			run_cmd += my::tostring(t_(i)) + ","; 
		}
		title   += my::tostring(t_.back()) + "),";
		run_cmd += " -d:Jp 1 -u:tmax 10 -d";

		rst_file_->title(title,'-'); 
		rst_file_->change_text_onclick("run command",run_cmd);

		rst_file_->figure(dir_+filename_+"-pstricks.png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(dir_+filename_+"-pstricks.pdf")+RST::scale("200"));
	}
}

void HoneycombFree::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="honeycomb-free";
	display_results();
}
/*}*/
