#include "KagomeDirac.hpp"

KagomeDirac::KagomeDirac(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc):
	System(ref,N,m,n,M,bc),
	Kagome<double>(2,3,6,"kagome-dirac")
{
	if(status_==1){
		init_fermionic();

		system_info_.text("KagomeDirac : 3 sites per unit cell, pi-flux per hexagon,");
		system_info_.text("no flux per triangle");
	}
}

/*{method needed for running*/
void KagomeDirac::compute_H(){
	double t(1.0);
	H_.set(n_,n_,0);
	Matrix<int> nb;
	unsigned int s(0);
	for(unsigned int i(0);i<Lx_;i++){
		for(unsigned int j(0);j<Ly_;j++){
			/*site 0*/
			s = spuc_*(i + j*Lx_);
			nb = get_neighbourg(s);
			/*0-1*/H_(s,nb(0,0)) = nb(0,1)*t;
			/*0-1*/H_(s,nb(2,0)) = nb(2,1)*t;

			/*site 1*/
			s++;
			nb = get_neighbourg(s);
			/*0-1*/H_(s,nb(1,0)) = nb(1,1)*t;
			/*0-1*/H_(s,nb(3,0)) = -nb(3,1)*t;

			/*site 2*/
			s++;
			nb = get_neighbourg(s);
			/*0-1*/H_(s,nb(0,0)) = -nb(0,1)*t;
			/*0-1*/H_(s,nb(2,0)) = nb(2,1)*t;

			/*site 3*/
			s++;
			nb = get_neighbourg(s);
			/*0-1*/H_(s,nb(0,0)) = nb(0,1)*t;
			/*0-1*/H_(s,nb(2,0)) = -nb(2,1)*t;

			/*site 4*/
			s++;
			nb = get_neighbourg(s);
			/*0-1*/H_(s,nb(1,0)) = nb(1,1)*t;
			/*0-1*/H_(s,nb(3,0)) = -nb(3,1)*t;

			/*site 5*/
			s++;
			nb = get_neighbourg(s);
			/*0-1*/H_(s,nb(0,0)) = nb(0,1)*t;
			/*0-1*/H_(s,nb(2,0)) = nb(2,1)*t;
		}
	}
	H_ += H_.transpose();
}

void KagomeDirac::create(){
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
void KagomeDirac::lattice(){
	Matrix<int> nb;
	double x0;
	double x1;
	double y0;
	double y1;
	double ll(1.0);
	double ex(4.0*ll);
	double exy(2.0*ll*cos(2.0*M_PI/6.0));
	double ey(2.0*ll*sin(2.0*M_PI/6.0));
	std::string color("black");

	PSTricks ps("./","lattice");
	ps.add("\\begin{pspicture}(-1,-1)(16,10)%"+filename_);
	Matrix<double> cell(4,2);
	cell(0,0) = 0.0;
	cell(0,1) = 0.0;
	cell(1,0) = ex;
	cell(1,1) = 0.0;
	cell(2,0) = ex + exy;
	cell(2,1) = ey;
	cell(3,0) = exy;
	cell(3,1) = ey;
	ps.polygon(cell,"linewidth=1pt,linecolor=red");
	cell(1,0)*=Lx_;
	cell(2,0) = Lx_*ex + Ly_*exy;
	cell(2,1)*=Ly_;
	cell(3,0)*=Ly_;
	cell(3,1)*=Ly_;
	ps.polygon(cell,"linewidth=1pt,linecolor=red,linestyle=dashed");

	unsigned int s;
	for(unsigned int i(0);i<Lx_;i++) {
		for(unsigned int j(0);j<Ly_;j++) {
			/*site 0*/
			s = spuc_*(i+j*Lx_);
			nb = get_neighbourg(s);
			x0 = 0.2+i*ex+j*exy;
			y0 = 0.1+j*ey; 
			ps.put(x0-0.2,y0+0.2,tostring(s));
			x1 = x0+ll;
			y1 = y0;
			if(H_(s,nb(0,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*0-1*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0-ll;
			if(H_(s,nb(2,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*0-1*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);

			/*site 1*/
			s++;
			nb = get_neighbourg(s);
			x0 = x0+ll;
			double x3(x0+ll);
			double y3(y0);
			ps.put(x0+0.2,y0+0.2,tostring(s));
			x1 = x0+ll*cos(4.0*M_PI/6.0);
			y1 = y0+ll*sin(4.0*M_PI/6.0);
			if(H_(s,nb(1,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*1-2*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(10.0*M_PI/6.0);
			y1 = y0+ll*sin(10.0*M_PI/6.0);
			if(H_(s,nb(3,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*1-2*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);

			/*site 2*/
			s++;
			nb = get_neighbourg(s);
			x0 = x0+ll*cos(4.0*M_PI/6.0);
			y0 = y0+ll*sin(4.0*M_PI/6.0);
			ps.put(x0+0.2,y0,tostring(s));
			x1 = x0+ll*cos(2.0*M_PI/6.0);
			y1 = y0+ll*sin(2.0*M_PI/6.0);
			if(H_(s,nb(0,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*2-0*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(8.0*M_PI/6.0);
			y1 = y0+ll*sin(8.0*M_PI/6.0);
			if(H_(s,nb(2,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*2-0*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);

			/*site 3*/
			s++;
			nb = get_neighbourg(s);
			x0 = x3;
			y0 = y3;
			ps.put(x0-0.2,y0+0.2,tostring(s));
			x1 = x0+ll;
			if(H_(s,nb(0,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*2-0*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0-ll;
			if(H_(s,nb(2,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*2-0*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);

			/*site 4*/
			s++;
			nb = get_neighbourg(s);
			x0 = x0+ll;
			ps.put(x0+0.2,y0+0.2,tostring(s));
			x1 = x0+ll*cos(4.0*M_PI/6.0);
			y1 = y0+ll*sin(4.0*M_PI/6.0);
			if(H_(s,nb(1,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*2-0*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(10.0*M_PI/6.0);
			y1 = y0+ll*sin(10.0*M_PI/6.0);
			if(H_(s,nb(3,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*2-0*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);

			/*site 4*/
			s++;
			nb = get_neighbourg(s);
			x0 = x0+ll*cos(4.0*M_PI/6.0);
			y0 = y0+ll*sin(4.0*M_PI/6.0);
			ps.put(x0+0.2,y0,tostring(s));
			x1 = x0+ll*cos(2.0*M_PI/6.0);
			y1 = y0+ll*sin(2.0*M_PI/6.0);
			if(H_(s,nb(0,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*2-0*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(8.0*M_PI/6.0);
			y1 = y0+ll*sin(8.0*M_PI/6.0);
			if(H_(s,nb(2,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*2-0*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
		}
	}
	ps.add("\\end{pspicture}");
	ps.save(true,true);
}

void KagomeDirac::check(){
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
	
	//BandStructure<double> bs(H_,Lx_,Ly_,spuc_,bc_);
	lattice();
}
/*}*/

/*{method needed for analysing*/
std::string KagomeDirac::extract_level_7(){
	rst_file_ = new RSTFile(info_+path_+dir_,filename_);
	unsigned int nruns;
	unsigned int tmax;

	(*read_)>>nruns>>tmax;
	/* the +1 is the averages over all runs */
	for(unsigned int i(0);i<nruns+1;i++){ 
		(*read_)>>E_>>corr_>>long_range_corr_;
	}
	(*jd_write_)("energy per site",E_);

	rst_file_->text(read_->get_header());
	rst_file_->save(false);
	delete rst_file_;
	rst_file_ = NULL;

	return filename_;
}

std::string KagomeDirac::extract_level_6(){
	unsigned int nof(0);
	(*read_)>>nof;
	save();
	for(unsigned int i(0);i<nof;i++){
		(*read_)>>E_;
		(*data_write_)<<M_(0)<<" "<<E_.get_x()<<" "<<E_.get_dx()<<" "<<ref_(0)<<ref_(1)<<ref_(2)<<IOFiles::endl;
	}
	(*jd_write_)("energy per site",E_);

	return filename_;
}

std::string KagomeDirac::extract_level_4(){
	(*read_)>>E_;
	(*jd_write_)("energy per site",E_);
	(*data_write_)<<M_(0)<<" "<<E_.get_x()<<" "<<E_.get_dx()<<" "<<ref_(0)<<ref_(1)<<ref_(2)<<IOFiles::endl;

	return filename_;
}
/*}*/
