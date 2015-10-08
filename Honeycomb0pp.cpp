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
					H_(i,nb(0,0))= nb(0,1)*th;
					H_(i,nb(1,0))= nb(1,1)*th;
					H_(i,nb(2,0))= nb(2,1)*td_;
				}break;
			case 2:
				{
					H_(i,nb(0,0))= nb(0,1)*td_;
					H_(i,nb(1,0))= nb(1,1)*th;
					H_(i,nb(2,0))= nb(2,1)*th;
				}break;
			case 4:
				{
					H_(i,nb(0,0))= nb(0,1)*th;
					H_(i,nb(1,0))= nb(1,1)*td_;
					H_(i,nb(2,0))= nb(2,1)*th;
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
	match(1) = 0;
	if(my::are_equal(x,match)){ return 1; }
	match(0) = 1.0/3.0;
	match(1) = 1.0/3.0;
	if(my::are_equal(x,match)){ return 2; }
	match(0) = 2.0/3.0;
	match(1) = 1.0/3.0;
	if(my::are_equal(x,match)){ return 3; }
	match(0) = 2.0/3.0;
	match(1) = 2.0/3.0;
	if(my::are_equal(x,match)){ return 4; }
	match(0) = 0;
	match(1) = 2.0/3.0;
	if(my::are_equal(x,match)){ return 5; }
	return 7;
}

Matrix<double> Honeycomb0pp::set_ab(){
	Matrix<double> tmp(2,2);
	tmp(0,0) = 1.0;
	tmp(1,0) = 1.0;
	tmp(0,1) = -1.0;
	tmp(1,1) = 2.0;
	return tmp;
}
/*}*/

/*{method needed for checking*/
void Honeycomb0pp::lattice(std::string const& path, std::string const& filename){
	compute_H();
	Matrix<double> e(2,2);
	e(0,0) = 1.0/3.0;
	e(1,0) = 1.0/3.0;
	e(0,1) = -1.0/2.0;
	e(1,1) = 1.0/2.0;
	Matrix<double> inv_e(2,2);
	inv_e(0,0) = e(1,1);
	inv_e(1,0) =-e(1,0);
	inv_e(0,1) =-e(0,1);
	inv_e(1,1) = e(0,0);
	inv_e/=(e(0,0)*e(1,1)-e(1,0)*e(0,1));

	Matrix<int> nb;
	std::string color("black");
	std::string linestyle("solid");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(path,filename);
	ps.begin(-4,-10,20,10,filename_);
	for(unsigned int i(0);i<n_;i+=2) {
		xy0 = get_pos_in_lattice(i);
		set_pos_LxLy(xy0);
		//set_in_basis(xy0);
		xy0 = (inv_e*LxLy_*xy0).chop();
		nb = get_neighbourg(i);

		if(nb(0,1)<0){
			color = "red";
			xy1 = xy0;
			xy1(0) += 0.5;
			xy1(1) -= 1.0;
			ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(nb(0,0)));
		} else {
			color = "black";
			xy1 = get_pos_in_lattice(nb(0,0));
			set_pos_LxLy(xy1);
			//set_in_basis(xy1);
			xy1 = inv_e*LxLy_*xy1;
		}
		xy1 = xy1.chop();
		if(H_(i,nb(0,0))>0){ linestyle = "solid"; }
		else { linestyle = "dashed"; }
		/*x-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth=1pt,linecolor="+color+",linestyle="+linestyle);

		if(nb(1,1)<0){
			color = "red";
			xy1 = xy0;
			xy1(0) -= 0.5;
			xy1(1) += 1.0;
			ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(nb(1,0)));
		} else {
			color = "black";
			xy1 = get_pos_in_lattice(nb(1,0));
			set_pos_LxLy(xy1);
			//set_in_basis(xy1);
			xy1 = inv_e*LxLy_*xy1;
		}
		xy1 = xy1.chop();
		if(H_(i,nb(1,0))>0){ linestyle = "solid"; }
		else { linestyle = "dashed"; }
		/*y-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth=1pt,linecolor="+color+",linestyle="+linestyle);

		if(nb(2,1)<0){
			color = "red";
			xy1 = xy0;
			xy1(0) -= 0.5;
			xy1(1) -= 1.0;
			ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(nb(2,0)));
		} else {
			color = "black";
			xy1 = get_pos_in_lattice(nb(2,0));
			set_pos_LxLy(xy1);
			//set_in_basis(xy1);
			xy1 = inv_e*LxLy_*xy1;
		}
		xy1 = xy1.chop();
		if(H_(i,nb(2,0))>0){ linestyle = "solid"; }
		else { linestyle = "dashed"; }
		/*y-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth=1pt,linecolor="+color+",linestyle="+linestyle);
	}

	for(unsigned int i(0);i<n_;i++) {
		xy0 = get_pos_in_lattice(i);
		set_pos_LxLy(xy0);
		//set_in_basis(xy0);
		xy0 = (LxLy_*xy0).chop();
		xy0 = inv_e*xy0;
		ps.put(xy0(0)-0.20,xy0(1)+0.15,my::tostring(i));
	}

	Vector<double> Lx(2);
	Lx(0) = LxLy_(0,0);
	Lx(1) = LxLy_(1,0);
	Lx = inv_e*Lx;
	Vector<double> Ly(2);
	Ly(0) = LxLy_(0,1);
	Ly(1) = LxLy_(1,1);
	Ly = inv_e*Ly;

	Matrix<double> polygon(4,2);
	polygon(0,0)=0;
	polygon(0,1)=0;
	polygon(1,0)=Lx(0);
	polygon(1,1)=Lx(1);
	polygon(2,0)=Lx(0)+Ly(0);
	polygon(2,1)=Lx(1)+Ly(1);
	polygon(3,0)=Ly(0);
	polygon(3,1)=Ly(1);
	for(unsigned int i(0);i<polygon.row();i++){ polygon(i,0) -= 1; }
	ps.polygon(polygon,"linecolor=green");

	Vector<double> a(2);
	a(0) = ab_(0,0);
	a(1) = ab_(1,0);
	a = inv_e*a;
	Vector<double> b(2);
	b(0) = ab_(0,1);
	b(1) = ab_(1,1);
	b = inv_e*b;

	polygon(0,0)=0;
	polygon(0,1)=0;
	polygon(1,0)=a(0);
	polygon(1,1)=a(1);
	polygon(2,0)=a(0)+b(0);
	polygon(2,1)=a(1)+b(1);
	polygon(3,0)=b(0);
	polygon(3,1)=b(1);
	for(unsigned int i(0);i<polygon.row();i++){ polygon(i,0) -= 1; }
	ps.polygon(polygon,"linecolor=blue");

	ps.end(true,true,true);
}

void Honeycomb0pp::check(){
	///*{debug 1*/
	//Matrix<int> nb;
	//for(unsigned int i(0);i<n_;i++){
	//nb = get_neighbourg(i);
	//std::cout<<i<<" ";
	//for(unsigned int j(0);j<z_;j++){
	//std::cout<<nb(j,0)<<" ";
	//}
	//std::cout<<std::endl;
	//}
	//std::cout<<"<<<<<<<<<<<<<<"<<std::endl;
	//nb = get_neighbourg(2);
	//std::cout<<nb<<std::endl;
	///*}*/
	///*{debug 2*/
	//double t(1.0);
	//Matrix<int> nb;
	//Matrix<double> Ttest(n_,n_,0);
	//for(unsigned int s(0);s<n_;s++){
	//nb = get_neighbourg(s);
	//for(unsigned int i(0);i<z_;i++){
	//Ttest(s,nb(i,0)) = t;
	//}
	//}
	//for(unsigned int i(0);i<n_;i++){
	//for(unsigned int j(0);j<n_;j++){
	//if(std::abs(Ttest(i,j)-std::abs(H_(i,j)))>0.2){
	//std::cout<<i<<" "<<j<<std::endl;
	//}
	//}
	//}
	///*}*/
	///*{debug 3*/
	//unsigned int k(0);
	//for(unsigned int i(0);i<n_;i++){
	//for(unsigned int j(0);j<n_;j++){
	//if(H_(i,j)!=0){
	//k++;
	//std::cout<<i<<" "<<j<<" "<<H_(i,j)<<std::endl;
	//}
	//}
	//}
	//std::cout<<k<<" "<<links_.row()<<std::endl;
	///*}*/
	///*{debug 4*/
	//Matrix<int> nb;
	//for(unsigned int s(0);s<n_;s++){
	//nb = get_neighbourg(s);
	//for(unsigned int i(0);i<z_;i++){
	//if(nb(i,1)<0){std::cout<<s<<" "<<nb(i,0)<<std::endl;}
	//}
	//}
	///*}*/

	std::cout<<"<<<<<<<<<<<<<<"<<std::endl;
	std::cout<<0<<" "<<get_site_in_ab(0)<<std::endl;
	std::cout<<1<<" "<<get_site_in_ab(1)<<std::endl;
	std::cout<<12<<" "<<get_site_in_ab(12)<<std::endl;
	std::cout<<9<<" "<<get_site_in_ab(9)<<std::endl;
	std::cout<<3<<" "<<get_site_in_ab(3)<<std::endl;
	lattice("./","lattice");
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
		(*read_)>>E_>>corr_>>lr_corr_;
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
