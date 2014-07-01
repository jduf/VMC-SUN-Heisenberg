#include "KagomeFermi.hpp"

KagomeFermi::KagomeFermi(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc):
	System(ref,N,m,n,M,bc),
	Kagome<double>(1,1,3,"kagome-fermi")
{
	init_fermionic();
	compute_T();
	
	rst_.text("KagomeFermi : All hopping term are identical, no flux, 3 sites per unit cell");
}

KagomeFermi::~KagomeFermi(){}

void KagomeFermi::compute_T(){
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
			/*0-1*/T_(s,nb(3,0)) = nb(3,1)*t;

			/*site 2*/
			s++;
			nb = get_neighbourg(s);
			/*0-1*/T_(s,nb(0,0)) = nb(0,1)*t;
			/*0-1*/T_(s,nb(2,0)) = nb(2,1)*t;
		}
	}
	T_ += T_.transpose();
}

void KagomeFermi::compute_P(Matrix<double>& Px, Matrix<double>& Py){
	std::cerr<<"KagomeFermi::compute_P : undefined method"<<Px<<Py<<std::endl;
	Px.set(n_,n_,0);
	Py.set(n_,n_,0);
	for(unsigned int j(0);j<Ly_;j++){
		for(unsigned int i(0);i<Lx_-1;i++){
			for(unsigned int k(0);k<3;k++){
				Px(3*i + 3*j*Lx_ + k, 3*i + 3*j*Lx_ + k+3) = 1.0;
			}
		}
		for(unsigned int k(0);k<3;k++){
			Px(3*(Lx_-1) + 3*j*Lx_ +k,3*j*Lx_ + k) = bc_;
		}
	}
	for(unsigned int i(0);i<Lx_;i++){
		for(unsigned int j(0);j<Ly_-1;j++){
			for(unsigned int k(0);k<3;k++){
				Py(3*i + 3*j*Lx_ + k, 3*i + 3*(j+1)*Lx_ + k) = 1.0;
			}
		}
		for(unsigned int k(0);k<3;k++){
			Py(3*i + 3*(Ly_-1)*Lx_ + k, 3*i + k) = bc_;
		}
	}
}

void KagomeFermi::create(){
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

void KagomeFermi::lattice(){
	Matrix<int> nb;
	double x0;
	double x1;
	double y0;
	double y1;
	double ll(1.0);
	double ex(2.0*ll);
	double exy(2.0*ll*cos(2.0*M_PI/6.0));
	double ey(2.0*ll*sin(2.0*M_PI/6.0));
	std::string color("black");

	Matrix<double> xy(4,2);
	xy(0,0) = 0.0;
	xy(0,1) = 0.0;
	xy(1,0) = ex;
	xy(1,1) = 0.0;
	xy(2,0) = ex + exy;
	xy(2,1) = ey;
	xy(3,0) = exy;
	xy(3,1) = ey;
	PSTricks ps("./","lattice");
	ps.add("\\begin{pspicture}(-1,-1)(16,10)%"+filename_);
	ps.polygon(xy,"linewidth=1pt,linecolor=red");
	unsigned int s;
	for(unsigned int i(0);i<Lx_;i++) {
		for(unsigned int j(0);j<Ly_;j++) {
			/*site 0*/
			s = spuc_*(i+j*Lx_);
			nb = get_neighbourg(s);
			x0 = 0.2+i*ex+j*exy;
			/*0.05 is there so there is no problem with latex and it shows
			 * better which sites are in the unit cell*/
			y0 = 0.1+j*ey; 
			ps.put(x0-0.2,y0+0.2,tostring(s));
			x1 = x0+ll;
			y1 = y0;
			if(T_(s,nb(0,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*0-1*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0-ll;
			if(T_(s,nb(2,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*0-1*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);

			/*site 1*/
			s++;
			nb = get_neighbourg(s);
			x0 = x0+ll;
			ps.put(x0+0.2,y0+0.2,tostring(s));
			x1 = x0+ll*cos(4.0*M_PI/6.0);
			y1 = y0+ll*sin(4.0*M_PI/6.0);
			if(T_(s,nb(1,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*1-2*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(10.0*M_PI/6.0);
			y1 = y0+ll*sin(10.0*M_PI/6.0);
			if(T_(s,nb(3,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*1-2*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);

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
			/*2-0*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(8.0*M_PI/6.0);
			y1 = y0+ll*sin(8.0*M_PI/6.0);
			if(T_(s,nb(2,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*2-0*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
		}
	}
	ps.add("\\end{pspicture}");
	ps.save(true,true);
}

void KagomeFermi::check(){
	//Matrix<double> Px;
	//Matrix<double> Py;
	//compute_P(Px,Py);
	//BandStructure<double> bs(T_,Px,Py);
	lattice();
}
