#include "KagomeDirac.hpp"

KagomeDirac::KagomeDirac(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc):
	System(ref,N,m,n,M,bc),
	Kagome<double>(1,1,6,"kagome-dirac")
{
	if(status_==1){
		init_fermionic();
		compute_T();

		rst_.text("KagomeDirac : 3 sites per unit cell, pi-flux per hexagon,");
		rst_.text("no flux per triangle");
	}
}

/*{method needed for running*/
void KagomeDirac::compute_T(){
	double t(1.0);
	T_.set(n_,n_,0);
	Matrix<int> nb;
	unsigned int s(0);
	for(unsigned int i(0);i<Lx_;i++){
		for(unsigned int j(0);j<Ly_;j++){
			/*site 0*/
			s = spuc_*(i + j*Lx_);
			nb = get_neighbourg(s);
			/*0-1*/T_(s,nb(0,0)) = nb(0,1)*t;
			/*0-1*/T_(s,nb(2,0)) = nb(2,1)*t;

			/*site 1*/
			s++;
			nb = get_neighbourg(s);
			/*0-1*/T_(s,nb(1,0)) = nb(1,1)*t;
			/*0-1*/T_(s,nb(3,0)) = -nb(3,1)*t;

			/*site 2*/
			s++;
			nb = get_neighbourg(s);
			/*0-1*/T_(s,nb(0,0)) = -nb(0,1)*t;
			/*0-1*/T_(s,nb(2,0)) = nb(2,1)*t;

			/*site 3*/
			s++;
			nb = get_neighbourg(s);
			/*0-1*/T_(s,nb(0,0)) = nb(0,1)*t;
			/*0-1*/T_(s,nb(2,0)) = -nb(2,1)*t;

			/*site 4*/
			s++;
			nb = get_neighbourg(s);
			/*0-1*/T_(s,nb(1,0)) = nb(1,1)*t;
			/*0-1*/T_(s,nb(3,0)) = -nb(3,1)*t;

			/*site 5*/
			s++;
			nb = get_neighbourg(s);
			/*0-1*/T_(s,nb(0,0)) = nb(0,1)*t;
			/*0-1*/T_(s,nb(2,0)) = nb(2,1)*t;
		}
	}
	T_ += T_.transpose();
}

void KagomeDirac::create(){
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);

	diagonalize_T();
	for(unsigned int c(0);c<N_;c++){
		if(!is_degenerate(c)){
			EVec_[c].set(n_,M_(c));
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = T_(i,j);
				}
			}
		}
	}
}
/*}*/

/*{method needed for checking*/
void KagomeDirac::compute_P(Matrix<double>& Px, Matrix<double>& Py){
	Px.set(n_,n_,0);
	Py.set(n_,n_,0);
	unsigned int tmp;
	for(unsigned int j(0);j<Ly_;j++){
		for(unsigned int i(0);i<Lx_-1;i++){
			tmp = spuc_*(i + j*Lx_);
			for(unsigned int k(0);k<spuc_;k++){
				Px(tmp+k, tmp+k+spuc_) = 1.0;
			}
		}
		tmp = spuc_*((Lx_-1) + j*Lx_);
		for(unsigned int k(0);k<spuc_;k++){
			Px(tmp+k,spuc_*j*Lx_ + k) = bc_;
		}
	}
	for(unsigned int i(0);i<Lx_;i++){
		for(unsigned int j(0);j<Ly_-1;j++){
			tmp = spuc_*(i + j*Lx_);
			for(unsigned int k(0);k<spuc_;k++){
				Py(tmp+k, tmp+spuc_*Lx_+k) = 1.0;
			}
		}
		tmp = spuc_*(i + (Ly_-1)*Lx_);
		for(unsigned int k(0);k<spuc_;k++){
			Py(tmp+k, spuc_*i+k) = bc_;
		}
	}
}

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
			if(T_(s,nb(0,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*0-1*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0-ll;
			if(T_(s,nb(2,0))>0){ color = "green"; }
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
			if(T_(s,nb(1,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*1-2*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(10.0*M_PI/6.0);
			y1 = y0+ll*sin(10.0*M_PI/6.0);
			if(T_(s,nb(3,0))>0){ color = "green"; }
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
			if(T_(s,nb(0,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*2-0*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(8.0*M_PI/6.0);
			y1 = y0+ll*sin(8.0*M_PI/6.0);
			if(T_(s,nb(2,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*2-0*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);

			/*site 3*/
			s++;
			nb = get_neighbourg(s);
			x0 = x3;
			y0 = y3;
			ps.put(x0-0.2,y0+0.2,tostring(s));
			x1 = x0+ll;
			if(T_(s,nb(0,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*2-0*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0-ll;
			if(T_(s,nb(2,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*2-0*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);

			/*site 4*/
			s++;
			nb = get_neighbourg(s);
			x0 = x0+ll;
			ps.put(x0+0.2,y0+0.2,tostring(s));
			x1 = x0+ll*cos(4.0*M_PI/6.0);
			y1 = y0+ll*sin(4.0*M_PI/6.0);
			if(T_(s,nb(1,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*2-0*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(10.0*M_PI/6.0);
			y1 = y0+ll*sin(10.0*M_PI/6.0);
			if(T_(s,nb(3,0))>0){ color = "green"; }
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
			if(T_(s,nb(0,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*2-0*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(8.0*M_PI/6.0);
			y1 = y0+ll*sin(8.0*M_PI/6.0);
			if(T_(s,nb(2,0))>0){ color = "green"; }
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
			//if(std::abs(Ttest(i,j)-std::abs(T_(i,j)))>0.2){
				//std::cout<<i<<" "<<j<<std::endl;
			//}
		//}
	//}
	///*}*/
	///*{debug 3*/
	//unsigned int k(0);
	//for(unsigned int i(0);i<n_;i++){
		//for(unsigned int j(0);j<n_;j++){
			//if(T_(i,j)!=0){
				//k++;
				//std::cout<<i<<" "<<j<<" "<<T_(i,j)<<std::endl;
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
	
	//Matrix<double> Px;
	//Matrix<double> Py;
	//compute_P(Px,Py);
	//BandStructure<double> bs(T_,Px,Py);
	lattice();
}
/*}*/
