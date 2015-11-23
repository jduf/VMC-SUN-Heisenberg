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
	double th(1.0);
	H_.set(n_,n_,0);
	Matrix<int> nb;
	unsigned int s(0);
	for(unsigned int i(0);i<n_;i+=2){
		s = get_site_in_ab(i);
		nb = get_neighbourg(i);
		switch(s){
			case 0:
				{
					H_(i,nb(0,0))= nb(0,1)*td_;
					H_(i,nb(1,0))= nb(1,1)*th;
					H_(i,nb(2,0))= nb(2,1)*th;
				}break;
			case 2:
				{
					H_(i,nb(0,0))= nb(0,1)*th;
					H_(i,nb(1,0))= nb(1,1)*td_;
					H_(i,nb(2,0))= nb(2,1)*th;
				}break;
			case 4:
				{
					H_(i,nb(0,0))= nb(0,1)*th;
					H_(i,nb(1,0))= nb(1,1)*th;
					H_(i,nb(2,0))= nb(2,1)*td_;
				}break;
			default:{ std::cerr<< __PRETTY_FUNCTION__<<" : undefined site in unit cell"<<std::endl; }break;
		}
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

void Honeycomb0pp::save_param(IOFiles& w) const{
	w.write("td/th (ratio of the hopping parameters)",td_);
}

unsigned int Honeycomb0pp::match_pos_in_ab(Vector<double> const& x) const{
	Vector<double> match(2,0);
	if(my::are_equal(x,match)){ return 0; }
	match(0) = 1.0/3.0;
	if(my::are_equal(x,match)){ return 1; }
	match(1) = 1.0/3.0;
	if(my::are_equal(x,match)){ return 2; }
	match(0) = 0;
	match(1) = 2.0/3.0;
	if(my::are_equal(x,match)){ return 3; }
	match(0) = 2.0/3.0;
	if(my::are_equal(x,match)){ return 4; }
	match(1) = 1.0/3.0;
	if(my::are_equal(x,match)){ return 5; }
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

	Matrix<int> nb;
	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	double t;
	PSTricks ps(info_+path_+dir_,filename_);
	ps.begin(-20,-10,20,10,filename_);

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

	for(unsigned int i(0);i<n_;i++){
		xy0 = get_pos_in_lattice(i);
		set_pos_LxLy(xy0);
		xy0 = (LxLy_*xy0).chop();
		ps.put(xy0(0)-0.20,xy0(1)+0.15,my::tostring(i));

		if(!(i%2)){
			nb = get_neighbourg(i);

			t = H_(i,nb(0,0));
			if(std::abs(t)>1e-5){
				xy1 = get_pos_in_lattice(nb(0,0));
				set_pos_LxLy(xy1);
				xy1 = LxLy_*xy1;
				if((xy0-xy1).norm_squared()>1.0001){
					linestyle = "dashed"; 
					xy1 = xy0;
					xy1(0) += dir_nn_(0,0);
					xy1(1) += dir_nn_(0,1);
					ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(nb(0,0)));
				} else { linestyle = "solid";  }
				xy1 = xy1.chop();

				if(t>0){ color = "blue";}
				else { color = "red"; }
				linewidth = my::tostring(std::abs(t))+"mm";
				/*(+x)-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
			}

			t = H_(i,nb(1,0));
			if(std::abs(t)>1e-5){
				xy1 = get_pos_in_lattice(nb(1,0));
				set_pos_LxLy(xy1);
				xy1 = LxLy_*xy1;
				if((xy0-xy1).norm_squared()>1.0001){
					linestyle = "dashed"; 
					xy1 = xy0;
					xy1(0) += dir_nn_(1,0);
					xy1(1) += dir_nn_(1,1);
					ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(nb(1,0)));
				} else { linestyle = "solid";  }
				xy1 = xy1.chop();

				if(t>0){ color = "blue";}
				else { color = "red"; }
				linewidth = my::tostring(std::abs(t))+"mm";
				/*(+y)-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
			}

			t = H_(i,nb(2,0));
			if(std::abs(t)>1e-5){
				xy1 = get_pos_in_lattice(nb(2,0));
				set_pos_LxLy(xy1);
				xy1 = LxLy_*xy1;
				if((xy0-xy1).norm_squared()>1.0001){
					linestyle = "dashed"; 
					xy1 = xy0;
					xy1(0) += dir_nn_(2,0);
					xy1(1) += dir_nn_(2,1);
					ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(nb(2,0)));
				} else { linestyle = "solid";  }
				xy1 = xy1.chop();

				if(t>0){ color = "blue";}
				else { color = "red"; }
				linewidth = my::tostring(std::abs(t))+"mm";
				/*(-y)-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
			}
		}
	}
	ps.end(true,true,true);
}

void Honeycomb0pp::check(){
	check_lattice();
	info_ ="";
	path_ ="";
	dir_ ="./";
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
