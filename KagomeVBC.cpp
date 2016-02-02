#include "KagomeVBC.hpp"

KagomeVBC::KagomeVBC(System const& s):
	System(s),
	Kagome<std::complex<double> >(set_ab(),9,"kagome-vbc")
{
	if(status_==2){
		init_fermionic();

		system_info_.text("KagomeVBC : 9 sites per unit cell, pi-flux through 1/3 of the hexagon");
		system_info_.text("and -pi/6-flux through all triangles, so the total flux is null");
	}
}

/*{method needed for running*/
void KagomeVBC::compute_H(){
	H_.set(n_,n_,0);
	double t(1.0);
	double phi(M_PI/6.0);

	unsigned int s0(0);
	unsigned int s1(0);
	unsigned int ab0(0);
	unsigned int ab1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		ab0 = obs_[0](i,5);
		ab1 = obs_[0](i,6);
		if((ab0==3 && ab1==5) || (ab0==4 && ab1==6) || (ab0==5 && ab1==2) || (ab0==7 && ab1==8) || (ab0==8 && ab1==0)){ H_(s0,s1) = std::polar((obs_[0](i,4)?bc_:1)*t,-phi); }
		else { H_(s0,s1) = std::polar((obs_[0](i,4)?bc_:1)*t,phi); }
	}
	H_ += H_.conjugate_transpose();
}

void KagomeVBC::create(){
	compute_H();
	diagonalize(true);
	if(status_==2){
		for(unsigned int c(0);c<N_;c++){
			EVec_[c].set(n_,M_(c));
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = H_(i,j);
				}
			}
		}
	}
	compute_H();
}

Matrix<double> KagomeVBC::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 3.0;
	tmp(1,0) =-sqrt(3.0);
	tmp(0,1) = 0.0;
	tmp(1,1) = 2*sqrt(3.0);
	return tmp;
}

unsigned int KagomeVBC::match_pos_in_ab(Vector<double> const& x) const {
	Vector<double> match(2,0);
	if(my::are_equal(x,match)){ return 0; }
	match(0) = 0.5;
	if(my::are_equal(x,match)){ return 1; }
	match(0) = 0.0;
	match(1) = 0.5;
	if(my::are_equal(x,match)){ return 2; }
	double a(1.0/3.0);
	double b(1.0/6.0);
	match(0) = a;
	match(1) = b; 
	if(my::are_equal(x,match)){ return 3; }
	match(0) += a;
	match(1) += b; 
	if(my::are_equal(x,match)){ return 4; }
	match(0) -= 0.5;
	if(my::are_equal(x,match)){ return 5; }
	match(0) += 2*a;
	match(1) += 2*b;
	if(my::are_equal(x,match)){ return 6; }
	match(0) -= 0.5;
	if(my::are_equal(x,match)){ return 7; }
	match(0) += a;
	match(1) += b;
	if(my::are_equal(x,match)){ return 8; }
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return 9;
}
/*}*/

/*{method needed for checking*/
void KagomeVBC::display_results(){
	compute_H();

	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	std::complex<double> t;
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
		if(std::abs(t)>1e-4){
			if((xy0-xy1).norm_squared()>1.0001){
				linestyle = "dashed";
				xy1 = xy0;
				switch(obs_[0](i,5)){
					case 2:
						{
							xy1(0) += dir_nn_(1,0);
							xy1(1) += dir_nn_(1,1);
						}break;
					case 4:
						{
							xy1(0) += dir_nn_(0,0);
							xy1(1) += dir_nn_(0,1);
						}break;
					case 6:
						{
							xy1(0) += dir_nn_(2,0);
							xy1(1) += dir_nn_(2,1);
						}break;
					case 7:
						{
							xy1(0) += dir_nn_(2,0);
							xy1(1) += dir_nn_(2,1);
						}break;
					case 8:
						{
							xy1(0) += dir_nn_(1,0);
							xy1(1) += dir_nn_(1,1);
						}break;
					default:
						{ std::cerr<<__PRETTY_FUNCTION__<<" : boundary not defined"<<std::endl; }
				}
				xy1 = xy1.chop();
				ps.put(xy1(0)-0.20,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
			} else { linestyle = "solid"; }

			if(t.imag()>0){ color = "blue"; }
			else          { color = "red"; }
			linewidth = my::tostring(std::abs(t))+"mm";
			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}
		ps.put(xy0(0)-0.20,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}");
	}
	ps.end(true,true,true);
}

void KagomeVBC::check(){
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
	///*}*/
	///*{debug 2*/
	//Matrix<int> nb;
	//double t(1.0);
	//Matrix<double> Ttest(n_,n_,0);
	//for(unsigned int s(0);s<n_;s++){
	//nb = get_neighbourg(s);
	//for(unsigned int i(0);i<z_;i++){ Ttest(s,nb(i,0)) = t; }
	//}
	//for(unsigned int i(0);i<n_;i++){
	//for(unsigned int j(0);j<n_;j++){
	//if(std::abs(Ttest(i,j)-my::norm_squared(H_(i,j)))>0.2){
	//std::cout<<i<<" "<<j<<std::endl;
	//}
	//}
	//}
	///*}*/
	///*{debug 3*/
	//unsigned int k(0);
	//for(unsigned int i(0);i<n_;i++){
	//for(unsigned int j(0);j<n_;j++){
	//if(my::norm_squared(H_(i,j))!=0){
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
	//if(nb(i,1)<0){ std::cout<<s<<" "<<nb(i,0)<<std::endl; }
	//}
	//}
	///*}*/

	//plot_band_structure();
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="kagome-vbc";
	display_results();
}
/*}*/

/*{method needed for analysing*/
std::string KagomeVBC::extract_level_7(){
	rst_file_ = new RSTFile(info_+path_+dir_,filename_);
	unsigned int nruns;
	unsigned int tmax;

	(*read_)>>nruns>>tmax;
	obs_.push_back(Observable(*read_));
	/* the +1 is the averages over all runs */
	jd_write_->write("E",obs_[0][0]);

	rst_file_->text(read_->get_header());
	rst_file_->save(false,true);
	delete rst_file_;
	rst_file_ = NULL;

	return filename_;
}

std::string KagomeVBC::extract_level_6(){
	unsigned int nof(0);
	(*read_)>>nof;
	save_input(*jd_write_);
	for(unsigned int i(0);i<nof;i++){
		(*read_)>>obs_[0][0];
		(*data_write_)<<M_(0)<<" "<<obs_[0][0]<<" "<<ref_(0)<<ref_(1)<<ref_(2)<<IOFiles::endl;
	}
	jd_write_->write("E",obs_[0][0]);

	return filename_;
}

std::string KagomeVBC::extract_level_4(){
	(*read_)>>obs_[0];
	jd_write_->write("E",obs_[0][0]);
	(*data_write_)<<M_(0)<<" "<<obs_[0][0]<<" "<<ref_(0)<<ref_(1)<<ref_(2)<<IOFiles::endl;

	return filename_;
}
/*}*/
