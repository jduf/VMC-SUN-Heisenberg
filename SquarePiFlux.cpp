#include "SquarePiFlux.hpp"

SquarePiFlux::SquarePiFlux(System const& s):
	System(s,3),
	Square<std::complex<double> >(set_ab(),(N_/m_==2?2:0),"square-csl")
{
	if(status_==2){
		init_fermionic();

		system_info_.text("Chiral spin liquid : pi-flux per plaquette");
	}
}

/*{method needed for running*/
void SquarePiFlux::compute_H(){
	double phi(M_PI/4.0);
	H_.set(n_,n_,0);
	Matrix<int> nb;
	unsigned int s(0);
	for(unsigned int i(0);i<n_;i++){
		s = get_site_in_ab(i);
		nb = get_neighbourg(i);
		switch(s){ 
			case 0:
				{
					H_(i,nb(0,0)) = std::polar(double(nb(0,1)),phi);
					H_(i,nb(1,0)) = std::polar(double(nb(1,1)),-phi);
				}break;
			case 1:
				{
					H_(i,nb(0,0)) = std::polar(double(nb(0,1)),-phi);
					H_(i,nb(1,0)) = std::polar(double(nb(1,1)),phi);
				}break;
			default:{ std::cerr<<"void SquarePiFlux::compute_H() : undefined site in unit cell"<<std::endl; }break;
		}
	}
	std::cerr<<"SquarePiFlux : compute_EVec : new use of polar, check that it is correct"<<std::endl;
	std::cerr<<"                            : modified the flux disposition..."<<std::endl;
	std::cerr<<"it seems that std::polar is not very stable for std::polar(1,-pi)=(0,1e-6)"<<std::endl;
	H_ += H_.trans_conj(); 
}

void SquarePiFlux::create(){
	compute_H();
	diagonalize(true);
	for(unsigned int c(0);c<N_;c++){
		EVec_[c].set(n_,M_(c));
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_(c);j++){
				EVec_[c](i,j) = H_(i,j);
			}
		}
	}
}

unsigned int SquarePiFlux::match_pos_in_ab(Vector<double> const& x) const{
	Vector<double> match(2,0);
	if(my::are_equal(x,match)){ return 0; }
	match(0) = 0.5;
	match(1) = 0;
	if(my::are_equal(x,match)){ return 1; }
	return 2;
}

Matrix<double> SquarePiFlux::set_ab(){
	Matrix<double> tmp(2,2);
	tmp(0,0) = 2;
	tmp(1,0) = 0;
	tmp(0,1) = 1;
	tmp(1,1) = 1;
	return tmp;
}
/*}*/

/*{method needed for checking*/
void SquarePiFlux::lattice(std::string const& path){
	compute_H();
	Matrix<int> nb;
	std::string color("black");
	std::string linestyle("solid");
	std::string arrow("-");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(path,"lattice");
	ps.add("\\begin{pspicture}(-9,-10)(16,10)%"+filename_);
	for(unsigned int i(0);i<n_;i++) {
		xy0 = get_pos_in_lattice(i);
		set_pos_LxLy(xy0);
		set_in_basis(xy0);
		xy0 = (LxLy_*xy0).chop();
		ps.put(xy0(0)-0.20,xy0(1)+0.15,my::tostring(i));
		nb = get_neighbourg(i);

		if(nb(0,1)<0){
			color = "red";
			xy1 = xy0;
			xy1(0) += 1.0;
			ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(nb(0,0)));
		} else {
			color = "black";
			if(i+1==xloop_){
				xy1 = xy0;
				xy1(0) += 1.0;
			} else {
				xy1 = get_pos_in_lattice(nb(0,0));
				set_pos_LxLy(xy1);
				set_in_basis(xy1);
				xy1 = (LxLy_*xy1).chop();
			}
		}
		if(H_(i,nb(0,0)).imag()>0){ arrow = "->"; }
		else { arrow = "<-"; }
		/*x-link*/ ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1), "linewidth=1pt,linecolor="+color+",linestyle="+linestyle);

		if(nb(1,1)<0){
			color = "red";
			xy1 = xy0;
			xy1(1) += 1.0;
			ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(nb(1,0)));
		} else {
			color = "black";
			xy1 = get_pos_in_lattice(nb(1,0));
			set_pos_LxLy(xy1);
			set_in_basis(xy1);
			xy1 = (LxLy_*xy1).chop();
		}
		if(H_(i,nb(1,0)).imag()>0){ arrow = "->"; }
		else { arrow = "<-"; }
		/*y-link*/ ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1), "linewidth=1pt,linecolor="+color+",linestyle="+linestyle);
	}

	Matrix<double> polygon(4,2);
	polygon(0,0)=0;
	polygon(0,1)=0;
	polygon(1,0)=LxLy_(0,0);
	polygon(1,1)=LxLy_(1,0);
	polygon(2,0)=LxLy_(0,0)+LxLy_(0,1);
	polygon(2,1)=LxLy_(1,0)+LxLy_(1,1);
	polygon(3,0)=LxLy_(0,1);
	polygon(3,1)=LxLy_(1,1);
	for(unsigned int i(0);i<polygon.row();i++){ polygon(i,0) -= 0.5; }
	ps.polygon(polygon,"linecolor=green");

	polygon(0,0)=0;
	polygon(0,1)=0;
	polygon(1,0)=ab_(0,0);
	polygon(1,1)=ab_(1,0);
	polygon(2,0)=ab_(0,0)+ab_(0,1);
	polygon(2,1)=ab_(1,0)+ab_(1,1);
	polygon(3,0)=ab_(0,1);
	polygon(3,1)=ab_(1,1);
	for(unsigned int i(0);i<polygon.row();i++){ 
		polygon(i,0) -= 0.2;
		polygon(i,1) -= 0.1;
	}
	ps.polygon(polygon,"linecolor=blue");

	ps.add("\\end{pspicture}");
	ps.save(true,true,true);
}

void SquarePiFlux::check(){
	//Matrix<int> nb;
	//for(unsigned int s(0);s<n_;s++){
	//nb = get_neighbourg(s);
	//for(unsigned int i(0);i<z_;i++){
	//std::cout<<s<<" "<<nb(i,0)<<" "<<nb(i,1)<<std::endl;
	//}
	//}
	lattice("./");
}
/*}*/

/*{method needed for analysing*/
std::string SquarePiFlux::extract_level_7(){
	rst_file_ = new RSTFile(info_+path_+dir_,filename_);

	unsigned int nruns;
	unsigned int tmax;

	(*read_)>>nruns>>tmax;
	(*data_write_)<<"% E dE 0|1"<<IOFiles::endl;
	/* the +1 is the averages over all runs */
	for(unsigned int i(0);i<nruns+1;i++){ 
		(*read_)>>E_>>corr_>>lr_corr_;
		(*data_write_)<<E_.get_x()<<" "<<E_.get_dx()<<" "<<(i<nruns?true:false)<<IOFiles::endl;
	}
	jd_write_->write("energy per site",E_);

	rst_file_->text(read_->get_header());
	rst_file_->save(false);
	delete rst_file_;
	rst_file_ = NULL;

	return filename_;
}

std::string SquarePiFlux::extract_level_3(){
	(*read_)>>E_;
	(*data_write_)<<n_<<" "<<E_<<" "<<bc_<<IOFiles::endl;

	return filename_;
}
/*}*/

//{//csl for Vishvanath (uses majorana representation)
//for(unsigned int i(0); i< Ly_; i++){
//for(unsigned int j(0); j< Lx_; j++){
//if(j+1 == Lx_){// x hopping
//H(i*Lx_ , i*Lx_ + j) = t;
//if(i % 2 == 0){
//T(i*Lx_ , i*Lx_ + j) = bc_*t;
//} else {
//T(i*Lx_ , i*Lx_ + j) = -bc_*t;
//}
//} else {
//H( i*Lx_ + j , i*Lx_ + j + 1) = t; 
//if(i % 2 == 0){
//T( i*Lx_ + j , i*Lx_ + j + 1) = t; 
//} else {
//T( i*Lx_ + j , i*Lx_ + j + 1) = -t; 
//}
//}
//if(i+1 == Ly_ ){// y hopping
//H(j, i*Lx_ + j) = t;
//T(j, i*Lx_ + j) = bc_*t;
//} else{
//H(i*Lx_ + j, (i+1)*Lx_ + j) = t;
//T(i*Lx_ + j, (i+1)*Lx_ + j) = t;
//}
//}
//}
//H += H.transpose();
//T += T.transpose();
//} 
