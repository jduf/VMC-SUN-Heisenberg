#include "Honeycomb0pp.hpp"

Honeycomb0pp::Honeycomb0pp(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc, double td):
	System(ref,N,m,n,M,bc),
	Honeycomb<double>(4,3,6,"honeycomb0pp"),
	td_(td)
{
	if(status_==1){
		init_fermionic();
		compute_T();

		filename_ += "-td" + tostring(td_);
		rst_.text("Honeycomb0pp : 6 sites per unit cell, in the center hexagon there is a 0-flux,");
		rst_.text("if td<0, the two other hexagons contain a pi-flux, if td>0, their flux is 0");
	}
}

/*{method needed for running*/
void Honeycomb0pp::compute_T(){
	double th(1.0);
	T_.set(n_,n_,0);
	Matrix<int> nb;
	unsigned int s(0);
	for(unsigned int i(0);i<Lx_;i++){
		for(unsigned int j(0);j<Ly_;j++){
			/*site 0*/
			s = spuc_*(i+j*Lx_);
			nb = get_neighbourg(s);
			T_(s,nb(0,0)) = nb(0,1)*th;
			T_(s,nb(1,0)) = nb(1,1)*td_;
			T_(s,nb(2,0)) = nb(2,1)*th;

			/*site 2*/
			s+=2;
			nb = get_neighbourg(s);
			T_(s,nb(0,0)) = nb(0,1)*th;
			T_(s,nb(1,0)) = nb(1,1)*td_;
			T_(s,nb(2,0)) = nb(2,1)*th;
			
			/*site 2*/
			s+=2;
			nb = get_neighbourg(s);
			T_(s,nb(0,0)) = nb(0,1)*th;
			T_(s,nb(1,0)) = nb(1,1)*td_;
			T_(s,nb(2,0)) = nb(2,1)*th;
		}
	}
	T_ += T_.transpose();
}

void Honeycomb0pp::create(){
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);

	diagonalize_T();
	for(unsigned int c(0);c<N_;c++){
		if(!is_degenerate(c)){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = T_(i,j);
				}
			}
		}
	}
}

void Honeycomb0pp::save(IOFiles& w) const{
	GenericSystem<double>::save(w);
	w("td/th (ratio of the hopping parameters)",td_);
}
/*}*/

/*{method needed for checking*/
void Honeycomb0pp::compute_P(Matrix<double>& Px, Matrix<double>& Py){
	Px.set(n_,n_,0);
	Py.set(n_,n_,0);
	for(unsigned int j(0);j<Ly_;j++){
		for(unsigned int i(0);i<Lx_-1;i++){
			for(unsigned int k(0);k<spuc_;k++){
				Px(spuc_*(i + j*Lx_) + k, spuc_*(i + j*Lx_)+ k+spuc_) = 1.0;
			}
		}
		for(unsigned int k(0);k<spuc_;k++){
			Px(spuc_*((Lx_-1) + j*Lx_) +k,spuc_*j*Lx_ + k) = bc_;
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

void Honeycomb0pp::lattice(){
	Matrix<int> nb;
	double angle_p(2.0*M_PI/6.0);
	double angle_n(-2.0*M_PI/6.0);
	double x0;
	double x1;
	double y0;
	double y1;
	double ll(1.0);
	double ex(2.0*ll*(1+cos(angle_p)));
	double exy(3.0*ll*cos(angle_p));
	double ey(3.0*ll*sin(angle_p));
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

	unsigned int s(0);
	for(unsigned int i(0);i<Lx_;i++){
		for(unsigned int j(0);j<Ly_;j++){
			/*site 0*/
			s = spuc_*(i+j*Lx_);
			nb = get_neighbourg(s);
			x0 = 0.5*ll+1.5*ll*cos(angle_p)+i*ex+j*exy;
			y0 = 1.5*ll*sin(angle_p)+j*ey; 
			ps.put(x0+0.2,y0,tostring(s));
			x1 = x0+ll*cos(angle_p);
			y1 = y0+ll*sin(angle_p);
			if(T_(s,nb(0,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*0-1*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0-ll;
			y1 = y0;
			if(T_(s,nb(1,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*0-3*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(angle_n);
			y1 = y0+ll*sin(angle_n);
			if(T_(s,nb(2,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*0-5*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);

			/*site 1*/
			s++;
			x0 = x0 + ll*cos(angle_p);
			y0 = y0 + ll*sin(angle_p);
			ps.put(x0+0.2,y0-0.2,tostring(s));

			/*site 2*/
			s++;
			nb = get_neighbourg(s);
			x0 = x0 + ll;
			ps.put(x0-0.2,y0-0.2,tostring(s));
			x1 = x0+ll*cos(angle_n);
			y1 = y0+ll*sin(angle_n);
			if(T_(s,nb(0,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*0-3*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(angle_p);
			y1 = y0+ll*sin(angle_p);
			if(T_(s,nb(1,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*0-5*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0-ll;
			y1 = y0;
			if(T_(s,nb(2,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*0-1*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);

			/*site 3*/
			s++;
			x0 = x0 + ll*cos(angle_n);
			y0 = y0 + ll*sin(angle_n);
			ps.put(x0-0.2,y0,tostring(s));

			/*site 4*/
			s++;
			nb = get_neighbourg(s);
			x0 = x0 - ll*cos(angle_p) ;
			y0 = y0 - ll*sin(angle_p);
			ps.put(x0-0.2,y0+0.2,tostring(s));
			x1 = x0-ll;
			y1 = y0;
			if(T_(s,nb(0,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*0-3*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(angle_n);
			y1 = y0+ll*sin(angle_n);
			if(T_(s,nb(1,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*0-5*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			x1 = x0+ll*cos(angle_p);
			y1 = y0+ll*sin(angle_p);
			if(T_(s,nb(2,0))>0){ color = "green"; }
			else { color = "blue"; }
			/*0-1*/	ps.line("-",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);

			/*site 5*/
			s++;
			x0 = x0 - ll;
			ps.put(x0+0.2,y0+0.2,tostring(s));
		}
	}
	ps.add("\\end{pspicture}");
	ps.save(true,true);
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
