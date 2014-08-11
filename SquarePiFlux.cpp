#include "SquarePiFlux.hpp"

SquarePiFlux::SquarePiFlux(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc):
	System(ref,N,m,n,M,bc),
	Square<std::complex<double> >(1,1,N,"square-csl")
{
	if(status_==1){
		init_fermionic();

		system_info_.text("Chiral spin liquid, with 2pi/N flux per plaquette");
		std::cout<<"check everything"<<std::endl;
	}
}

/*{method needed for running*/
void SquarePiFlux::compute_H(){
	double t(1.0);
	double phi(2.0*M_PI/spuc_);
	H_.set(n_,n_,0);
	Matrix<int> nb;
	unsigned int s(0);
	for(unsigned int i(0);i<Ly_;i++){
		for(unsigned int j(0);j<Lx_;j++){
			s = spuc_*(i + j*Lx_);
			for(unsigned int k(0);k<spuc_;k++){
				nb = get_neighbourg(s);
				H_(s,nb(0,0)) = nb(0,1)*t;
				H_(s,nb(1,0)) = std::polar(nb(1,1)*t,k*phi);;
				s++;
			}
		}
	}
	std::cerr<<"SquarePiFlux : compute_EVec : new use of polar, check that it is correct"<<std::endl;
	std::cerr<<"                            : modified the flux disposition..."<<std::endl;
	std::cerr<<"it seems that std::polar is not very stable for std::polar(1,-pi)=(0,1e-6)"<<std::endl;
	H_ += H_.trans_conj(); 
}

void SquarePiFlux::create(){
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);

	compute_H();
	diagonalize_H(H_);
	for(unsigned int c(0);c<N_;c++){
		EVec_[c].set(n_,M_(c));
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_(c);j++){
				EVec_[c](i,j) = H_(i,j);
			}
		}
	}
}
/*}*/

/*{method needed for checking*/
void SquarePiFlux::lattice(){
	Matrix<int> nb;
	double x0;
	double x1;
	double y0;
	double y1;
	double ll(1.0);
	double ex(spuc_*ll);
	double ey(ll);
	std::string color;

	PSTricks ps("./","lattice");
	ps.add("\\begin{pspicture}(-1,-1)(16,10)%"+filename_);
	unsigned int s;
	for(unsigned int i(0);i<Lx_;i++) {
		for(unsigned int j(0);j<Ly_;j++) {
			s = spuc_*(i+j*Lx_);
			x0 = i*ex;
			y0 = j*ey;
			for(unsigned int k(0);k<spuc_;k++){
				nb = get_neighbourg(s);
				ps.put(x0-0.2,y0+0.2,tostring(s));
				x1 = x0+ll;
				y1 = y0;
				if(real(H_(s,nb(0,0)))>0){ color = "green"; }
				else { color = "blue"; }
				/*x-link*/ ps.line("-",x0,y0,x1,y1, "linewidth=1pt,linecolor="+color);
				x1 = x0;
				y1 = y0+ll;
				if(real(H_(s,nb(1,0)))>0){ color = "green"; }
				else { color = "blue"; }
				/*y-link*/ ps.line("-",x0,y0,x1,y1, "linewidth=1pt,linecolor="+color);
				x0 += ll;
				s++;
			}
		}
	}

	ps.frame(-0.5,-0.5,Lx_*ex-0.5,Ly_*ey-0.5,"linecolor=red");
	ps.frame(-0.5,-0.5,ex-0.5,ey-0.5,"linecolor=red,linestyle=dashed");
	ps.add("\\end{pspicture}");
	ps.save(true,true);
}

void SquarePiFlux::check(){
	Matrix<int> nb;
	for(unsigned int s(0);s<n_;s++){
		nb = get_neighbourg(s);
		for(unsigned int i(0);i<z_;i++){
			std::cout<<s<<" "<<nb(i,0)<<" "<<nb(i,1)<<std::endl;
		}
	}
	compute_H();
	lattice();
}
/*}*/

/*{method needed for analysing*/
std::string SquarePiFlux::extract_level_7(){
	rst_file_ = new RSTFile(info_+path_+dir_,filename_);

	unsigned int nruns;
	unsigned int tmax;

	(*read_)>>nruns>>tmax;
	data_write_->precision(10);
	(*data_write_)<<"% E dE 0|1"<<IOFiles::endl;
	/* the +1 is the averages over all runs */
	for(unsigned int i(0);i<nruns+1;i++){ 
		(*read_)>>E_>>corr_>>lr_corr_;
		(*data_write_)<<E_.get_x()<<" "<<E_.get_dx()<<" "<<(i<nruns?true:false)<<IOFiles::endl;
	}
	(*jd_write_)("energy per site",E_);

	rst_file_->text(read_->get_header());
	rst_file_->save(false);
	delete rst_file_;
	rst_file_ = NULL;

	return filename_;
}

std::string SquarePiFlux::extract_level_3(){
	(*read_)>>E_;
	(*data_write_)<<n_<<" "<<E_<<IOFiles::endl;

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
