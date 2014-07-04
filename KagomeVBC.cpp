#include "KagomeVBC.hpp"

KagomeVBC::KagomeVBC(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc):
	System(ref,N,m,n,M,bc),
	Kagome<std::complex<double> >(1,1,9,"kagome-vbc")
{
	if(status_==1){
		init_fermionic();
		compute_T();

		rst_.text("KagomeVBC : 9 sites per unit cell, pi-flux through 1/3 of the hexagon");
		rst_.text("and -pi/6-flux through all triangles, so the total flux is null");
	}
}

/*{method needed for running*/
void KagomeVBC::compute_T(){
	double t(1.0);
	double phi(M_PI/6.0);
	T_.set(n_,n_,0);
	Matrix<int> nb;
	unsigned int s(0);
	for(unsigned int j(0);j<Ly_;j++){
		for(unsigned int i(0);i<Lx_;i++){
			/*site 0*/
			s = spuc_*(i + j*Lx_);
			nb = get_neighbourg(s);
			/*0-1*/ T_(s, nb(0,0)) = std::polar(nb(0,1)*t,phi);
			/*0-8*/ T_(s, nb(1,0)) = std::polar(nb(1,1)*t,phi);
			/*0-6*/ T_(s, nb(2,0)) = std::polar(nb(2,1)*t,-phi);

			/*site 1*/
			s++;
			nb = get_neighbourg(s);
			/*1-2*/ T_(s, nb(0,0)) = std::polar(nb(0,1)*t,-phi);
			/*1-7*/ T_(s, nb(1,0)) = std::polar(nb(1,1)*t,phi); 
			/*1-8*/ T_(s, nb(2,0)) = std::polar(nb(2,1)*t,-phi); 

			/*site 2*/
			s++;
			nb = get_neighbourg(s);
			/*2-3*/ T_(s, nb(0,0)) = std::polar(nb(0,1)*t,phi);
			/*2-6*/ T_(s, nb(1,0)) = std::polar(nb(1,1)*t,phi);
			/*2-7*/ T_(s, nb(2,0)) = std::polar(nb(2,1)*t,-phi);

			/*site 3*/
			s++;
			nb = get_neighbourg(s);
			/*3-4*/ T_(s, nb(0,0)) = std::polar(nb(0,1)*t,-phi);
			/*3-8*/ T_(s, nb(1,0)) = std::polar(nb(1,1)*t,-phi);
			/*3-6*/ T_(s, nb(2,0)) = std::polar(nb(2,1)*t,-phi);

			/*site 4*/
			s++;
			nb = get_neighbourg(s);
			/*4-5*/ T_(s, nb(0,0)) = std::polar(nb(0,1)*t,phi);
			/*4-7*/ T_(s, nb(1,0)) = std::polar(nb(1,1)*t,-phi);
			/*4-8*/ T_(s, nb(2,0)) = std::polar(nb(2,1)*t,-phi);

			/*site 5*/
			s++;
			nb = get_neighbourg(s);
			/*5-0*/ T_(s, nb(0,0)) = std::polar(nb(0,1)*t,-phi);
			/*5-6*/ T_(s, nb(1,0)) = std::polar(nb(1,1)*t,-phi);
			/*5-7*/ T_(s, nb(2,0)) = std::polar(nb(2,1)*t,-phi);

		}
	}
	T_ += T_.trans_conj();
}

void KagomeVBC::create(){
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

void KagomeVBC::lattice(){
	Matrix<int> nb;
	double x0;
	double x1;
	double y0;
	double y1;
	double ll(1.0);
	double ex(4.0*ll*cos(M_PI/6.0));
	double exy(2.0*ll*cos(M_PI/6.0));
	double ey(3.0);
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
	cell(1,0)*= Lx_;
	cell(2,0) = Lx_*ex + Ly_*exy;
	cell(2,1)*= Ly_;
	cell(3,0)*= Ly_;
	cell(3,1)*= Ly_;
	ps.polygon(cell,"linewidth=1pt,linecolor=red,linestyle=dashed");

	unsigned int s;
	for(unsigned int i(0);i<Lx_;i++) {
		for(unsigned int j(0);j<Ly_;j++) {
			/*site 0*/
			s = spuc_*(i+j*Lx_);
			nb = get_neighbourg(s);
			x0 = ll*(0.5+sin(M_PI/6.0)) + i*ex + j*exy;
			/*0.05 is there so there is no problem with latex and it shows
			 * better which sites are in the unit cell*/
			y0 = 0.05 + ll/2.0 + j*ey; 
			ps.put(x0+0.2,y0+0.2,tostring(s));
			x1 = x0+ll*cos(3.0*M_PI/6.0);
			y1 = y0+ll*sin(3.0*M_PI/6.0);
			if(imag(T_(s,nb(0,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*0-1*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(9.0*M_PI/6.0);
			y1 = y0+ll*sin(9.0*M_PI/6.0);
			if(imag(T_(s,nb(2,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*0-6*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);


			/*site 1*/
			s++;
			nb = get_neighbourg(s);
			x0 = x0+ll*cos(3.0*M_PI/6.0);
			y0 = y0+ll*sin(3.0*M_PI/6.0);
			ps.put(x0+0.2,y0-0.2,tostring(s));
			x1 = x0+ll*cos(1.0*M_PI/6.0);
			y1 = y0+ll*sin(1.0*M_PI/6.0);
			if(imag(T_(s,nb(0,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*1-2*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(7.0*M_PI/6.0);
			y1 = y0+ll*sin(7.0*M_PI/6.0);
			double x8(x1);
			double y8(y1);
			if(imag(T_(s,nb(2,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*1-8*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);

			/*site 2*/
			s++;
			nb = get_neighbourg(s);
			x0 = x0+ll*cos(1.0*M_PI/6.0);
			y0 = y0+ll*sin(1.0*M_PI/6.0);
			ps.put(x0,y0-0.2,tostring(s));
			x1 = x0+ll*cos(-1.0*M_PI/6.0);
			y1 = y0+ll*sin(-1.0*M_PI/6.0);
			if(imag(T_(s,nb(0,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*2-3*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(5.0*M_PI/6.0);
			y1 = y0+ll*sin(5.0*M_PI/6.0);
			double x7(x1);
			double y7(y1);
			if(imag(T_(s,nb(2,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*2-7*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);

			/*site 3*/
			s++;
			nb = get_neighbourg(s);
			x0 = x0+ll*cos(-1.0*M_PI/6.0);
			y0 = y0+ll*sin(-1.0*M_PI/6.0);
			ps.put(x0-0.2,y0-0.2,tostring(s));
			x1 = x0+ll*cos(9.0*M_PI/6.0);
			y1 = y0+ll*sin(9.0*M_PI/6.0);
			if(imag(T_(s,nb(0,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*3-4*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(3.0*M_PI/6.0);
			y1 = y0+ll*sin(3.0*M_PI/6.0);
			double x6(x1);
			double y6(y1);
			if(imag(T_(s,nb(2,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*3-6*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);

			/*site 4*/
			s++;
			nb = get_neighbourg(s);
			x0 = x0+ll*cos(-3.0*M_PI/6.0);
			y0 = y0+ll*sin(-3.0*M_PI/6.0);
			ps.put(x0-0.2,y0+0.2,tostring(s));
			x1 = x0+ll*cos(7.0*M_PI/6.0);
			y1 = y0+ll*sin(7.0*M_PI/6.0);
			if(imag(T_(s,nb(0,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*4-5*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(1.0*M_PI/6.0);
			y1 = y0+ll*sin(1.0*M_PI/6.0);
			if(imag(T_(s,nb(2,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*4-8*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);

			/*site 5*/
			s++;
			nb = get_neighbourg(s);
			x0 = x0-ll*cos(1.0*M_PI/6.0);
			y0 = y0-ll*sin(1.0*M_PI/6.0);
			ps.put(x0,y0+0.2,tostring(s));
			x1 = x0+ll*cos(5.0*M_PI/6.0);
			y1 = y0+ll*sin(5.0*M_PI/6.0);
			if(imag(T_(s,nb(0,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*5-0*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(11.0*M_PI/6.0);
			y1 = y0+ll*sin(11.0*M_PI/6.0);
			if(imag(T_(s,nb(2,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*5-7*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);

			/*site 6*/
			s++;
			nb = get_neighbourg(s);
			x0 = x6;
			y0 = y6;
			ps.put(x0+0.2,y0-0.2,tostring(s));
			x1 = x0+ll*cos(1.0*M_PI/6.0);
			y1 = y0+ll*sin(1.0*M_PI/6.0);
			if(imag(T_(s,nb(0,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*6-5*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(7.0*M_PI/6.0);
			y1 = y0+ll*sin(7.0*M_PI/6.0);
			if(imag(T_(s,nb(2,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*6-2*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);

			/*site 7*/
			s++;
			nb = get_neighbourg(s);
			x0 = x7+ex;
			y0 = y7;
			ps.put(x0-0.2,y0-0.2,tostring(s));
			x1 = x0+ll*cos(3.0*M_PI/6.0);
			y1 = y0+ll*sin(3.0*M_PI/6.0);
			if(imag(T_(s,nb(0,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*7-1*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(9.0*M_PI/6.0);
			y1 = y0+ll*sin(9.0*M_PI/6.0);
			if(imag(T_(s,nb(2,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*7-4*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);

			/*site 8*/
			s++;
			nb = get_neighbourg(s);
			x0 = x8+ex;
			y0 = y8;
			ps.put(x0,y0+0.2,tostring(s));
			x1 = x0+ll*cos(11.0*M_PI/6.0);
			y1 = y0+ll*sin(11.0*M_PI/6.0);
			if(imag(T_(s,nb(0,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*8-0*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(5.0*M_PI/6.0);
			y1 = y0+ll*sin(5.0*M_PI/6.0);
			if(imag(T_(s,nb(2,0)))>0){ color = "green";}
			else { color = "blue"; }
			/*8-3*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
		}
	}
	ps.add("\\end{pspicture}");
	ps.save(true,true);
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
			//if(std::abs(Ttest(i,j)-norm_squared(T_(i,j)))>0.2){
				//std::cout<<i<<" "<<j<<std::endl;
			//}
		//}
	//}
	///*}*/
	///*{debug 3*/
	//unsigned int k(0);
	//for(unsigned int i(0);i<n_;i++){
		//for(unsigned int j(0);j<n_;j++){
			//if(norm_squared(T_(i,j))!=0){
				//k++;
				//std::cout<<i<<" "<<j<<" "<<T_(i,j)<<std::endl;
			//}
		//}
	//}
	//std::cout<<k<<" "<<links_.row()<<std::endl;
	///*}*/
	/*{debug 4*/
	Matrix<int> nb;
	for(unsigned int s(0);s<n_;s++){
		nb = get_neighbourg(s);
		for(unsigned int i(0);i<z_;i++){
			if(nb(i,1)<0){std::cout<<s<<" "<<nb(i,0)<<std::endl;}
		}
	}
	/*}*/
	
	//Matrix<std::complex<double> > Px;
	//Matrix<std::complex<double> > Py;
	//compute_P(Px,Py);
	//BandStructure<std::complex<double> > bs(T_,Px,Py);
	lattice();
}
/*}*/
