#include "Honeycomb0pp.hpp"

Honeycomb0pp::Honeycomb0pp(System const& s, double td):
	System(s),
	Honeycomb<double>(set_ab(),6,"honeycomb0pp"),
	td_(td)
{
	if(status_==2){
		init_lattice();
		init_fermionic();

		filename_ += "-td" + my::tostring(td_);
		system_info_.text("Honeycomb0pp : 6 sites per unit cell, in the center hexagon there is a 0-flux,");
		system_info_.text("if td<0, the two other hexagons contain a pi-flux, if td>0, their flux is 0");
		system_info_.text("th is set to 1");
	}
}

/*{method needed for running*/
void Honeycomb0pp::compute_H(){
	H_.set(n_,n_,0);

	unsigned int s0(0);
	unsigned int s1(0);
	unsigned int ab0(0);
	unsigned int ab1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		ab0 = obs_[0](i,5);
		ab1 = obs_[0](i,6);
		if((ab0==0 && ab1==1) || (ab0==2 && ab1==3) || (ab0==4 && ab1==5)){ H_(s0,s1) = (obs_[0](i,4)?bc_:1)*td_; } 
		else { H_(s0,s1) = (obs_[0](i,4)?bc_:1); }
	}
	H_ += H_.transpose();
}

void Honeycomb0pp::create(){
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

void Honeycomb0pp::save_param(IOFiles& w) const {
	w.write("td/th (ratio of the hopping parameters)",td_);
}

Matrix<double> Honeycomb0pp::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 3.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 1.5;
	tmp(1,1) = 1.5*sqrt(3.0);
	return tmp;
}

unsigned int Honeycomb0pp::match_pos_in_ab(Vector<double> const& x) const {
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
/*}*/

/*{method needed for checking*/
void Honeycomb0pp::lattice(){
	compute_H();

	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(info_+path_+dir_,filename_);
	ps.begin(-20,-20,20,20,filename_);
	ps.polygon(lattice_corners_,"linecolor=green");

	double x_shift(-ab_(0,1)*3.0/2.0);
	double y_shift((-ab_(1,0)-ab_(1,1))/2);
	Matrix<double> polygon(4,2);
	polygon(0,0)=x_shift;
	polygon(0,1)=y_shift;
	polygon(1,0)=x_shift+ab_(0,0);
	polygon(1,1)=y_shift+ab_(1,0);
	polygon(2,0)=x_shift+ab_(0,0)+ab_(0,1);
	polygon(2,1)=y_shift+ab_(1,0)+ab_(1,1);
	polygon(3,0)=x_shift+ab_(0,1);
	polygon(3,1)=y_shift+ab_(1,1);
	ps.polygon(polygon,"linecolor=black");

	double t;
	double corr;
	unsigned int s0;
	unsigned int s1;
	std::string str;
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		xy0 = x_[s0];

		s1 = obs_[0](i,1);
		xy1 = x_[s1];

		//if(!(my::in_polygon(polygon.row(),polygon.ptr(),polygon.ptr()+polygon.row(),xy0(0),xy0(1)) || my::in_polygon(polygon.row(),polygon.ptr(),polygon.ptr()+polygon.row(),xy1(0),xy1(1))) ){
		t = H_(s0,s1);
		if(std::abs(t)>1e-4){
			if((xy0-xy1).norm_squared()>1.0001){
				linestyle = "dashed";
				xy1 = (xy0+dir_nn_[obs_[0](i,3)]).chop();
				ps.put(xy1(0)-0.20,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
			} else { 
				linestyle = "solid";  
				if(s0<s1){
					ps.put(xy0(0)-0.20,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}"); 
					ps.put(xy1(0)-0.20,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}"); 
				}
			}

			if(t>0){ color = "blue";}
			else   { color = "red"; }
			linewidth = my::tostring(std::abs(t))+"mm";
			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}
		//} else {
		if(obs_.size()>1){/*bound energy*/
			corr = obs_[1][obs_[0](i,2)].get_x();
			if(std::abs(corr)>1e-4){
				if(corr>0){ color = "blue"; }
				else      { color = "red"; }
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

void Honeycomb0pp::display_results(){
	lattice();
	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string run_cmd("./mc -s:wf ladder-free");
		run_cmd += " -u:N " + my::tostring(N_);
		run_cmd += " -u:m " + my::tostring(m_);
		run_cmd += " -u:n " + my::tostring(n_);
		run_cmd += " -i:bc "+ my::tostring(bc_);
		run_cmd += " -d:td " + my::tostring(td_); 
		run_cmd += " -d:Jp 1 -u:tmax 10 -d";
		rst_file_->change_text_onclick("run command",run_cmd);

		rst_file_->figure(dir_+filename_+".png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(dir_+filename_+".pdf")+RST::scale("200"));
		if(obs_[0].nval()){
			rst_file_->figure(relative_path+filename_+"-lr.png","long range correlations",RST::target(relative_path+filename_+"-lr.gp")+RST::scale("200"));
		} 
		if(obs_.size()==5){
			rst_file_->figure(relative_path+filename_+"-as.png","(anti)symmetric correlations",RST::target(relative_path+filename_+"-as.gp")+RST::scale("200"));
		}
	}
}

void Honeycomb0pp::check(){
	//obs_[0].print();
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="honeycomb-0pp";
	display_results();
	//std::cout<<get_site_in_ab(18)<<std::endl;
}
/*}*/

/*{method needed for analysing*/
std::string Honeycomb0pp::extract_level_7(){
	(*data_write_)<<td_<<" "<<obs_[0][0]<<IOFiles::endl;

	save_param(*jd_write_);
	save_output(*jd_write_);

	rst_file_ = new RSTFile(info_+path_+dir_,filename_);
	rst_file_->text(read_->get_header());
	rst_file_->save(false,true);
	delete rst_file_;
	rst_file_ = NULL;

	return my::tostring(td_);
}

std::string Honeycomb0pp::extract_level_6(){
	double td(0);
	unsigned int nof(0);
	(*read_)>>nof;
	int nobs;
	std::vector<Observable> obs;
	for(unsigned int i(0);i<nof;i++){
		(*read_)>>td>>nobs;
		obs.clear();
		for(int i(0);i<nobs;i++){ obs.push_back(Observable(*read_)); }
		if(obs[0][0].get_x()<get_energy().get_x()){
			obs_ = obs ;
			td_ = td;
		}
	}

	jd_write_->add_header()->nl();
	save_param(*jd_write_);
	save(*jd_write_);

	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp+="set xlabel '$\\frac{t_d}{t_h}$' offset 17,1.4";
	gp+="set y2label '$\\dfrac{E}{n}$'";
	gp+="plot '"+filename_+".dat' u 1:2:3 w e notitle";
	gp.save_file();
	gp.create_image(true,true);

	return filename_;
}
/*}*/
