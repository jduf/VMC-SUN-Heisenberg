#include "Honeycomb0pp.hpp"

Honeycomb0pp::Honeycomb0pp(System const& s, double td):
	System(s),
	Honeycomb<double>(set_ab(),6,"honeycomb0pp"),
	td_(td)
{
	if(status_==2){
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
	double th(1.0);

	unsigned int s0(0);
	unsigned int s1(0);
	unsigned int ab0(0);
	unsigned int ab1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		ab0 = get_site_in_ab(s0);
		ab1 = get_site_in_ab(s1);
		if((ab0==0 && ab1==1) || (ab0==2 && ab1==3) || (ab0==4 && ab1==5)){
			H_(s0,s1) = obs_[0](i,4)*td_; 
		} else { H_(s0,s1) = obs_[0](i,4)*th; }
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

Matrix<double> Honeycomb0pp::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 3.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 1.5;
	tmp(1,1) = 1.5*sqrt(3.0);
	return tmp;
}
/*}*/

/*{method needed for checking*/
void Honeycomb0pp::display_results(){
	compute_H();

	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	double t;
	PSTricks ps(info_+path_+dir_,filename_);
	ps.begin(-2,-20,40,20,filename_);

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

	unsigned int s0;
	unsigned int s1;
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		xy0 = get_pos_in_lattice(s0);

		s1 = obs_[0](i,1);
		xy1 = get_pos_in_lattice(s1);

		t = H_(s0,s1);
		if(std::abs(t)>1e-5){
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
					ps.put(xy0(0)-0.20,xy0(1)+0.15,my::tostring(s0)); 
					ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(s1)); 
				}
			}

			if(t>0){ color = "blue";}
			else { color = "red"; }
			linewidth = my::tostring(std::abs(t))+"mm";
			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}
	}
	ps.end(true,true,true);
}

void Honeycomb0pp::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="honeycomb0pp";
	display_results();
}
/*}*/

/*{method needed for analysing*/
std::string Honeycomb0pp::extract_level_7(){
	rst_file_ = new RSTFile(info_+path_+dir_,filename_);

	unsigned int nruns;
	unsigned int tmax;

	(*read_)>>nruns>>tmax;
	(*data_write_)<<"% td E dE 0|1"<<IOFiles::endl;
	/* the +1 is the averages over all runs */
	for(unsigned int i(0);i<nruns+1;i++){
		(*read_)>>E_>>obs_[0]>>obs_[1];
		(*data_write_)<<td_<<" "<<E_.get_x()<<" "<<E_.get_dx()<<" "<<(i<nruns?true:false)<<IOFiles::endl;
	}

	jd_write_->write("td/th (ratio of the hopping parameters)",td_);
	jd_write_->write("E",E_);

	rst_file_->text(read_->get_header());
	rst_file_->save(false,true);
	delete rst_file_;
	rst_file_ = NULL;

	return my::tostring(td_);
}

std::string Honeycomb0pp::extract_level_6(){
	Data<double> tmp_E;
	E_.set_x(1e33);
	double tmp_td(0);
	unsigned int nof(0);
	(*read_)>>nof;
	for(unsigned int i(0);i<nof;i++){
		(*read_)>>tmp_td>>tmp_E;
		if(tmp_E.get_x()<E_.get_x()){
			E_ = tmp_E;
			td_ = tmp_td;
		}
	}

	jd_write_->add_header()->nl();
	save_input(*jd_write_);
	jd_write_->write("energy per site",E_);

	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp+="set xlabel '$\\frac{t_d}{t_h}$' offset 17,1.4";
	gp+="set y2label '$\\dfrac{E}{n}$'";
	gp+="plot '"+filename_+".dat' u 1:($4==1?$2:1/0):3 w e t 'Independant measures',\\";
	gp+="     '"+filename_+".dat' u 1:($4==0?$2:1/0):3 w e t 'Mean'";
	gp.save_file();
	gp.create_image(true,true);

	return filename_;
}
/*}*/
