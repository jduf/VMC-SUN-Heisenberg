#include "SquarePiFlux.hpp"

SquarePiFlux::SquarePiFlux(System const& s):
	System(s),
	Square<std::complex<double> >((N_/m_==2?2:0),2,1,"square-csl")
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
			default:{ std::cerr<<__PRETTY_FUNCTION__<<" : undefined site in unit cell"<<std::endl; }break;
		}
	}
	std::cerr<<__PRETTY_FUNCTION__<<" : new use of polar, check that it is correct"<<std::endl;
	std::cerr<<__PRETTY_FUNCTION__<<" : modified the flux disposition..."<<std::endl;
	std::cerr<<__PRETTY_FUNCTION__<<" : it seems that std::polar is not very stable for std::polar(1,-pi)=(0,1e-6)"<<std::endl;
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
/*}*/

/*{method needed for checking*/
void SquarePiFlux::lattice(std::string const& path, std::string const& filename){
	compute_H();
	std::string color("black");
	std::string linestyle("solid");
	std::string arrow("-");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(path,filename);
	ps.begin(-9,-10,16,10,filename_);
	std::complex<double> t;
	unsigned int s0;
	unsigned int s1;
	double y_shift(4);
	for(unsigned int i(0);i<obs_[0].size();i++){
		s0 = obs_[0](i,0);
		xy0 = get_pos_in_lattice(s0);
		set_pos_LxLy(xy0);
		xy0 = (LxLy_*xy0).chop();

		s1 = obs_[0](i,1);
		xy1 = get_pos_in_lattice(s1);
		set_pos_LxLy(xy1);
		xy1 = (LxLy_*xy1).chop();

		if((xy0-xy1).norm_squared()<1.1){ linestyle = "solid"; }
		else {
			linestyle = "dashed";
			if(i%2 && xy1(1)<xy0(1)){
				xy1(0) = xy0(0);
				xy1(1) = xy0(1)+1.0;
			}
			if(!(i%2) && xy1(0)<xy0(0)){
				xy1(0) = xy0(0)+1.0;
				xy1(1) = xy0(1);
			}
			ps.put(xy1(0)-0.20,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
		}

		t = H_(s0,s1);
		if(i%2){
			ps.put(xy0(0)-0.20,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}");
		}

		if(std::abs(t)>1e-4){
			if(t.real()<0){ color = "red"; }
			else { color = "blue"; }

			if(t.imag()>0){ arrow = "->"; }
			else { arrow = "<-"; }
			xy0 = xy0.chop();
			xy1 = xy1.chop();
			ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1), "linewidth=1pt,linecolor="+color+",linestyle="+linestyle);
		}
	}

	double lr_corr;
	double rescale(obs_[1].size()?0.75/obs_[1][0].get_x():0);
	for(unsigned int i(0);i<obs_[1].size();i++){
		lr_corr = obs_[1][i].get_x()*rescale;
		if(std::abs(lr_corr)>1e-4){
			xy1 = get_pos_in_lattice(i);
			set_pos_LxLy(xy1);
			xy1 = (LxLy_*xy1).chop();
			xy1(1) -= 2*y_shift;

			if(i){
				if(lr_corr<0){ color = "red"; }
				else { color = "blue"; }
			} else { color = "black"; }

			ps.circle(xy1,std::abs(lr_corr),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
		}
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
	ps.polygon(polygon,"linecolor=black");
	ps.end(true,true,true);
}

void SquarePiFlux::check(){
	lattice("./","lattice");
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
		(*read_)>>E_>>obs_[0]>>obs_[1];
		(*data_write_)<<E_.get_x()<<" "<<E_.get_dx()<<" "<<(i<nruns?true:false)<<IOFiles::endl;
	}
	jd_write_->write("energy per site",E_);

	rst_file_->text(read_->get_header());
	rst_file_->save(false,true);
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
