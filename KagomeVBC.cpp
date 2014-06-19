#include "KagomeVBC.hpp"

KagomeVBC::KagomeVBC(unsigned int const& N, unsigned int const& n, unsigned int const& m, int const& bc, Vector<unsigned int> const& ref):
	System(N,n,m,bc,ref),
	Kagome<std::complex<double> >(1,1,9,"kagome-vbc")
{
	rst_.text("KagomeVBC : All hopping term are identical, therefore the unit cell contains only 3 sites");
}

KagomeVBC::~KagomeVBC(){}

void KagomeVBC::compute_T(){
	double t(1.0);
	double phi(M_PI/6.0);
	T_.set(n_,n_,0);
	Matrix<int> nb;
	/*0-flux per unit cell*/
	for(unsigned int j(0);j<Ly_;j++){
		for(unsigned int i(0);i<Lx_;i++){
			/*site 0*/
			nb = get_neighbourg(spuc_*(i + j*Lx_) + 0);
			/*0-1*/ T_(spuc_*(i + j*Lx_) + 0, nb(0,0)) = std::polar(nb(0,1)*t,phi);
			/*0-8*/ T_(spuc_*(i + j*Lx_) + 0, nb(1,0)) = std::polar(nb(1,1)*t,phi);
			/*0-6*/ T_(spuc_*(i + j*Lx_) + 0, nb(2,0)) = std::polar(nb(2,1)*t,-phi);

			/*site 1*/
			nb = get_neighbourg(spuc_*(i + j*Lx_) + 1);
			/*1-2*/ T_(spuc_*(i + j*Lx_) + 1, nb(0,0)) = std::polar(nb(0,1)*t,-phi);
			/*1-7*/ T_(spuc_*(i + j*Lx_) + 1, nb(1,0)) = std::polar(nb(1,1)*t,phi); 
			/*1-8*/ T_(spuc_*(i + j*Lx_) + 1, nb(2,0)) = std::polar(nb(2,1)*t,-phi); 

			/*site 2*/
			nb = get_neighbourg(spuc_*(i + j*Lx_) + 2);
			/*2-3*/ T_(spuc_*(i + j*Lx_) + 2, nb(0,0)) = std::polar(nb(0,1)*t,phi);
			/*2-6*/ T_(spuc_*(i + j*Lx_) + 2, nb(1,0)) = std::polar(nb(1,1)*t,phi);
			/*2-7*/ T_(spuc_*(i + j*Lx_) + 2, nb(2,0)) = std::polar(nb(2,1)*t,-phi);

			/*site 3*/
			nb = get_neighbourg(spuc_*(i + j*Lx_) + 3);
			/*3-4*/ T_(spuc_*(i + j*Lx_) + 3, nb(0,0)) = std::polar(nb(0,1)*t,-phi);
			/*3-8*/ T_(spuc_*(i + j*Lx_) + 3, nb(1,0)) = std::polar(nb(1,1)*t,-phi);
			/*3-6*/ T_(spuc_*(i + j*Lx_) + 3, nb(2,0)) = std::polar(nb(2,1)*t,-phi);

			/*site 4*/
			nb = get_neighbourg(spuc_*(i + j*Lx_) + 4);
			/*4-5*/ T_(spuc_*(i + j*Lx_) + 4, nb(0,0)) = std::polar(nb(0,1)*t,phi);
			/*4-7*/ T_(spuc_*(i + j*Lx_) + 4, nb(1,0)) = std::polar(nb(1,1)*t,-phi);
			/*4-8*/ T_(spuc_*(i + j*Lx_) + 4, nb(2,0)) = std::polar(nb(2,1)*t,-phi);

			/*site 5*/
			nb = get_neighbourg(spuc_*(i + j*Lx_) + 5);
			/*5-0*/ T_(spuc_*(i + j*Lx_) + 5, nb(0,0)) = std::polar(nb(0,1)*t,-phi);
			/*5-6*/ T_(spuc_*(i + j*Lx_) + 5, nb(1,0)) = std::polar(nb(1,1)*t,-phi);
			/*5-7*/ T_(spuc_*(i + j*Lx_) + 5, nb(2,0)) = std::polar(nb(2,1)*t,-phi);

		}
	}
	T_ += T_.trans_conj();
}

void KagomeVBC::compute_P(Matrix<std::complex<double> >& Px, Matrix<std::complex<double> >& Py){
	Px.set(n_,n_,0);
	Py.set(n_,n_,0);
	for(unsigned int j(0);j<Ly_;j++){
		for(unsigned int i(0);i<Lx_-1;i++){
			for(unsigned int k(0);k<spuc_;k++){
				Px(spuc_*(i + j*Lx_) + k, spuc_*(i + j*Lx_) + k + spuc_) = 1.0;
			}
		}
		for(unsigned int k(0);k<spuc_;k++){
			Px(spuc_*(Lx_-1 + j*Lx_) + k,spuc_*j*Lx_ + k) = bc_;
		}
	}
	for(unsigned int i(0);i<Lx_;i++){
		for(unsigned int j(0);j<Ly_-1;j++){
			for(unsigned int k(0);k<spuc_;k++){
				Py(spuc_*(i + j*Lx_) + k, spuc_*(i + (j+1)*Lx_) + k) = 1.0;
			}
		}
		for(unsigned int k(0);k<spuc_;k++){
			Py(spuc_*(i + (Ly_-1)*Lx_) + k, spuc_*i + k) = bc_;
		}
	}
}

void KagomeVBC::create(double const& x, unsigned int const& type){
	if(type!=2){std::cerr<<"KagomeVBC::create(double x, unsigned int const& type) : type unknown"<<x<<std::endl;}
	compute_T();
	diagonalize_T('S');
	for(unsigned int spin(0);spin<N_;spin++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_;j++){
				EVec_(i+spin*n_,j) = T_(i,j);
			}
		}
	}
}

void KagomeVBC::lattice(){
	PSTricks ps("./","lattice");
	Matrix<int> nb;
	double x0;
	double x1;
	double y0;
	double y1;
	double ex(4.0*cos(M_PI/6.0));
	double exy(2.0*cos(M_PI/6.0));
	double ey(3);
	std::string color("black");
	double ll(1.0);
	ps.add("\\begin{pspicture}(-1,-1)(16,10)%"+filename_);

	Matrix<double> xy(4,2);
	xy(0,0) = 0.0;
	xy(0,1) = 0.0;
	xy(1,0) = ex;
	xy(1,1) = 0.0;
	xy(2,0) = ex + exy;
	xy(2,1) = ey;
	xy(3,0) = exy;
	xy(3,1) = ey;
	ps.polygon(xy,"linewidth="+tostring(1)+"pt,linecolor=red");
	unsigned int s;

	for(unsigned int i(0);i<Lx_;i++) {
		for(unsigned int j(0);j<Ly_;j++) {
			/*site 0*/
			s = spuc_*(i+j*Lx_);
			nb = get_neighbourg(s);
			x0 = ll*(0.5+sin(M_PI/6.0)) + i*4.0*cos(M_PI/6.0) + j*2.0*cos(M_PI/6.0);
			/*0.05 is there so there is no problem with latex and it shows
			 * better which sites are in the unit cell*/
			y0 = 0.05 + ll/2.0 + j*3.0; 
			ps.put(x0+0.2,y0+0.2,tostring(s));
			x1 = x0+ll*cos(3.0*M_PI/6.0);
			y1 = y0+ll*sin(3.0*M_PI/6.0);
			if(imag(T_(s,nb(0,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*0-1*/	ps.line("->",x0,y0,x1,y1,"linewidth="+tostring(1)+"pt,linecolor="+color);
			x1 = x0+ll*cos(9.0*M_PI/6.0);
			y1 = y0+ll*sin(9.0*M_PI/6.0);
			if(imag(T_(s,nb(2,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*0-6*/	ps.line("->",x0,y0,x1,y1,"linewidth="+tostring(1)+"pt,linecolor="+color);


			/*site 1*/
			s = spuc_*(i+j*Lx_)+1;
			nb = get_neighbourg(s);
			x0 = x0+ll*cos(3.0*M_PI/6.0);
			y0 = y0+ll*sin(3.0*M_PI/6.0);
			ps.put(x0+0.2,y0-0.2,tostring(s));
			x1 = x0+ll*cos(1.0*M_PI/6.0);
			y1 = y0+ll*sin(1.0*M_PI/6.0);
			if(imag(T_(s,nb(0,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*1-2*/	ps.line("->",x0,y0,x1,y1,"linewidth="+tostring(1)+"pt,linecolor="+color);
			x1 = x0+ll*cos(7.0*M_PI/6.0);
			y1 = y0+ll*sin(7.0*M_PI/6.0);
			double x8(x1);
			double y8(y1);
			if(imag(T_(s,nb(2,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*1-8*/	ps.line("->",x0,y0,x1,y1,"linewidth="+tostring(1)+"pt,linecolor="+color);

			/*site 2*/
			s = spuc_*(i+j*Lx_)+2;
			nb = get_neighbourg(s);
			x0 = x0+ll*cos(1.0*M_PI/6.0);
			y0 = y0+ll*sin(1.0*M_PI/6.0);
			ps.put(x0,y0-0.2,tostring(spuc_*(i+j*Lx_)+2));
			x1 = x0+ll*cos(-1.0*M_PI/6.0);
			y1 = y0+ll*sin(-1.0*M_PI/6.0);
			if(imag(T_(s,nb(0,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*2-3*/	ps.line("->",x0,y0,x1,y1,"linewidth="+tostring(1)+"pt,linecolor="+color);
			x1 = x0+ll*cos(5.0*M_PI/6.0);
			y1 = y0+ll*sin(5.0*M_PI/6.0);
			double x7(x1);
			double y7(y1);
			if(imag(T_(s,nb(2,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*2-7*/	ps.line("->",x0,y0,x1,y1,"linewidth="+tostring(1)+"pt,linecolor="+color);

			/*site 3*/
			s = spuc_*(i+j*Lx_)+3;
			nb = get_neighbourg(s);
			x0 = x0+ll*cos(-1.0*M_PI/6.0);
			y0 = y0+ll*sin(-1.0*M_PI/6.0);
			ps.put(x0-0.2,y0-0.2,tostring(s));
			x1 = x0+ll*cos(9.0*M_PI/6.0);
			y1 = y0+ll*sin(9.0*M_PI/6.0);
			if(imag(T_(s,nb(0,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*3-4*/	ps.line("->",x0,y0,x1,y1,"linewidth="+tostring(1)+"pt,linecolor="+color);
			x1 = x0+ll*cos(3.0*M_PI/6.0);
			y1 = y0+ll*sin(3.0*M_PI/6.0);
			double x6(x1);
			double y6(y1);
			if(imag(T_(s,nb(2,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*3-6*/	ps.line("->",x0,y0,x1,y1,"linewidth="+tostring(1)+"pt,linecolor="+color);

			/*site 4*/
			s = spuc_*(i+j*Lx_)+4;
			nb = get_neighbourg(s);
			x0 = x0+ll*cos(-3.0*M_PI/6.0);
			y0 = y0+ll*sin(-3.0*M_PI/6.0);
			ps.put(x0-0.2,y0+0.2,tostring(s));
			x1 = x0+ll*cos(7.0*M_PI/6.0);
			y1 = y0+ll*sin(7.0*M_PI/6.0);
			if(imag(T_(s,nb(0,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*4-5*/	ps.line("->",x0,y0,x1,y1,"linewidth="+tostring(1)+"pt,linecolor="+color);
			x1 = x0+ll*cos(1.0*M_PI/6.0);
			y1 = y0+ll*sin(1.0*M_PI/6.0);
			if(imag(T_(s,nb(2,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*4-8*/	ps.line("->",x0,y0,x1,y1,"linewidth="+tostring(1)+"pt,linecolor="+color);

			/*site 5*/
			s = spuc_*(i+j*Lx_)+5;
			nb = get_neighbourg(s);
			x0 = x0-ll*cos(1.0*M_PI/6.0);
			y0 = y0-ll*sin(1.0*M_PI/6.0);
			ps.put(x0,y0+0.2,tostring(s));
			x1 = x0+ll*cos(5.0*M_PI/6.0);
			y1 = y0+ll*sin(5.0*M_PI/6.0);
			if(imag(T_(s,nb(0,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*5-0*/	ps.line("->",x0,y0,x1,y1,"linewidth="+tostring(1)+"pt,linecolor="+color);
			x1 = x0+ll*cos(11.0*M_PI/6.0);
			y1 = y0+ll*sin(11.0*M_PI/6.0);
			if(imag(T_(s,nb(2,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*5-7*/	ps.line("->",x0,y0,x1,y1,"linewidth="+tostring(1)+"pt,linecolor="+color);

			/*site 6*/
			s = spuc_*(i+j*Lx_)+6;
			nb = get_neighbourg(s);
			x0 = x6;
			y0 = y6;
			ps.put(x0+0.2,y0-0.2,tostring(s));
			x1 = x0+ll*cos(1.0*M_PI/6.0);
			y1 = y0+ll*sin(1.0*M_PI/6.0);
			if(imag(T_(s,nb(0,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*6-5*/	ps.line("->",x0,y0,x1,y1,"linewidth="+tostring(1)+"pt,linecolor="+color);
			x1 = x0+ll*cos(7.0*M_PI/6.0);
			y1 = y0+ll*sin(7.0*M_PI/6.0);
			if(imag(T_(s,nb(2,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*6-2*/	ps.line("->",x0,y0,x1,y1,"linewidth="+tostring(1)+"pt,linecolor="+color);

			/*site 7*/
			s = spuc_*(i+j*Lx_)+7;
			nb = get_neighbourg(s);
			x0 = x7+ex;
			y0 = y7;
			ps.put(x0-0.2,y0-0.2,tostring(s));
			x1 = x0+ll*cos(3.0*M_PI/6.0);
			y1 = y0+ll*sin(3.0*M_PI/6.0);
			if(imag(T_(s,nb(0,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*7-1*/	ps.line("->",x0,y0,x1,y1,"linewidth="+tostring(1)+"pt,linecolor="+color);
			x1 = x0+ll*cos(9.0*M_PI/6.0);
			y1 = y0+ll*sin(9.0*M_PI/6.0);
			if(imag(T_(s,nb(2,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*7-4*/	ps.line("->",x0,y0,x1,y1,"linewidth="+tostring(1)+"pt,linecolor="+color);

			/*site 8*/
			s = spuc_*(i+j*Lx_)+8;
			nb = get_neighbourg(s);
			x0 = x8+ex;
			y0 = y8;
			ps.put(x0,y0+0.2,tostring(s));
			x1 = x0+ll*cos(11.0*M_PI/6.0);
			y1 = y0+ll*sin(11.0*M_PI/6.0);
			if(imag(T_(s,nb(0,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*8-0*/	ps.line("->",x0,y0,x1,y1,"linewidth="+tostring(1)+"pt,linecolor="+color);
			x1 = x0+ll*cos(5.0*M_PI/6.0);
			y1 = y0+ll*sin(5.0*M_PI/6.0);
			if(imag(T_(s,nb(2,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*8-3*/	ps.line("->",x0,y0,x1,y1,"linewidth="+tostring(1)+"pt,linecolor="+color);
		}
	}
	ps.add("\\end{pspicture}");
	ps.save(true,true);
}

void KagomeVBC::check(){
	compute_T();
	Matrix<std::complex<double> > Px;
	Matrix<std::complex<double> > Py;
	compute_P(Px,Py);
	BandStructure<std::complex<double> > bs(T_,Px,Py);
}
