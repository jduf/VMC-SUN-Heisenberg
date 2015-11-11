#include "HoneycombSU4.hpp"

HoneycombSU4::HoneycombSU4(System const& s):
	System(s),
	Honeycomb<double>(set_ab(),6,"honeycomb-pi-flux")
{
	if(status_==2 && N_==4){
		init_fermionic();

		system_info_.text("SU(4) Honeycomb lattice");
		system_info_.item("pi-flux configuration");
		system_info_.item("4 sites per unit cell");
		system_info_.item("Lx x Ly unit cells");
	} else {
		std::cerr<<"HoneycombSU4 : the cluster is not a square"<<std::endl;
		std::cerr<<"HoneycombSU4 : or N!=4"<<std::endl;
	}
}

/*{method needed for running*/
void HoneycombSU4::compute_H(){
	H_.set(n_,n_,0);
	double th(1.0);
	double td(-1.0);
	//unsigned int i(0);
	/*old way to compute H*/
	//for(unsigned int l(0);l<Ly_;l++){
		//for(unsigned int c(0);c<Lx_;c++){
			//H_(i,i+1) = td; //0
			//if(l+1<Ly_){ H_(i,i+1+Lx_*4) = th; }
			//else { H_(i+1-l*Lx_*4,i) = th*bc_; }
			//if(c==0){ H_(i,i+Lx_*4-1) = th*bc_; }
			//else { H_(i-1,i) = th; }
			//i+=2;//2
			//H_(i,i-1) = th;
			//H_(i,i+1) = th; 
			//if(l==0){ H_(i,i+1+(Ly_-1)*Lx_*4) = th*bc_; }
			//else { H_(i,i+1-Lx_*4) = th; }
			//i+=2;//4
		//}
	//}
	/*new way to compute H (copy-pasted from Honeycomb0pp)*/
	Matrix<int> nb;
	unsigned int s(0);
	for(unsigned int i(0);i<n_;i+=2){
		s = get_site_in_ab(i);
		nb = get_neighbourg(i);
		switch(s){
			case 0:
				{
					H_(i,nb(0,0))= nb(0,1)*th;
					H_(i,nb(1,0))= nb(1,1)*th;
					H_(i,nb(2,0))= nb(2,1)*td;
				}break;
			case 2:
				{
					H_(i,nb(0,0))= nb(0,1)*td;
					H_(i,nb(1,0))= nb(1,1)*th;
					H_(i,nb(2,0))= nb(2,1)*th;
				}break;
			case 4:
				{
					H_(i,nb(0,0))= nb(0,1)*th;
					H_(i,nb(1,0))= nb(1,1)*td;
					H_(i,nb(2,0))= nb(2,1)*th;
				}break;
			default:{ std::cerr<< __PRETTY_FUNCTION__<<" : undefined site in unit cell"<<std::endl; }break;
		}
	}
	H_ += H_.transpose();
}

void HoneycombSU4::create(){
	compute_H();
	diagonalize(true);
	for(unsigned int c(0);c<N_;c++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_(c);j++){
				EVec_[c](i,j) = H_(i,j);
			}
		}
	}
}

/*copy-pasted from Honeycomb0pp*/
unsigned int HoneycombSU4::match_pos_in_ab(Vector<double> const& x) const{
	Vector<double> match(2,0);
	if(my::are_equal(x,match)){ return 0; }
	match(0) = 1.0/3.0;
	match(1) = 0;
	if(my::are_equal(x,match)){ return 1; }
	match(0) = 1.0/3.0;
	match(1) = 1.0/3.0;
	if(my::are_equal(x,match)){ return 2; }
	match(0) = 2.0/3.0;
	match(1) = 1.0/3.0;
	if(my::are_equal(x,match)){ return 3; }
	match(0) = 2.0/3.0;
	match(1) = 2.0/3.0;
	if(my::are_equal(x,match)){ return 4; }
	match(0) = 0;
	match(1) = 2.0/3.0;
	if(my::are_equal(x,match)){ return 5; }
	return 7;
}

/*copy-pasted from Honeycomb0pp*/
Matrix<double> HoneycombSU4::set_ab(){
	Matrix<double> tmp(2,2);
	tmp(0,0) = 1.0;
	tmp(1,0) = 1.0;
	tmp(0,1) = -1.0;
	tmp(1,1) = 2.0;
	return tmp;
}
/*}*/

/*{method needed for checking*/
void HoneycombSU4::check(){
	compute_H();
	std::cout<<H_.chop(1e-6)<<std::endl;
}

/*copy-pasted from Honeycomb0pp*/
void HoneycombSU4::display_results(){
	compute_H();
	Matrix<double> e(2,2);
	e(0,0) = 1.0/3.0;
	e(1,0) = 1.0/3.0;
	e(0,1) = -1.0/2.0;
	e(1,1) = 1.0/2.0;
	Matrix<double> inv_e(2,2);
	inv_e(0,0) = e(1,1);
	inv_e(1,0) =-e(1,0);
	inv_e(0,1) =-e(0,1);
	inv_e(1,1) = e(0,0);
	inv_e/=(e(0,0)*e(1,1)-e(1,0)*e(0,1));

	Matrix<int> nb;
	std::string color("black");
	std::string linestyle("solid");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(info_+path_+dir_,filename_);
	ps.begin(-4,-10,20,10,filename_);
	for(unsigned int i(0);i<n_;i+=2) {
		xy0 = get_pos_in_lattice(i);
		set_pos_LxLy(xy0);
		//set_in_basis(xy0);
		xy0 = (inv_e*LxLy_*xy0).chop();
		nb = get_neighbourg(i);

		if(nb(0,1)<0){
			color = "red";
			xy1 = xy0;
			xy1(0) += 0.5;
			xy1(1) -= 1.0;
			ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(nb(0,0)));
		} else {
			color = "black";
			xy1 = get_pos_in_lattice(nb(0,0));
			set_pos_LxLy(xy1);
			//set_in_basis(xy1);
			xy1 = inv_e*LxLy_*xy1;
		}
		xy1 = xy1.chop();
		if(H_(i,nb(0,0))>0){ linestyle = "solid"; }
		else { linestyle = "dashed"; }
		/*x-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth=1pt,linecolor="+color+",linestyle="+linestyle);

		if(nb(1,1)<0){
			color = "red";
			xy1 = xy0;
			xy1(0) -= 0.5;
			xy1(1) += 1.0;
			ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(nb(1,0)));
		} else {
			color = "black";
			xy1 = get_pos_in_lattice(nb(1,0));
			set_pos_LxLy(xy1);
			//set_in_basis(xy1);
			xy1 = inv_e*LxLy_*xy1;
		}
		xy1 = xy1.chop();
		if(H_(i,nb(1,0))>0){ linestyle = "solid"; }
		else { linestyle = "dashed"; }
		/*y-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth=1pt,linecolor="+color+",linestyle="+linestyle);

		if(nb(2,1)<0){
			color = "red";
			xy1 = xy0;
			xy1(0) -= 0.5;
			xy1(1) -= 1.0;
			ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(nb(2,0)));
		} else {
			color = "black";
			xy1 = get_pos_in_lattice(nb(2,0));
			set_pos_LxLy(xy1);
			//set_in_basis(xy1);
			xy1 = inv_e*LxLy_*xy1;
		}
		xy1 = xy1.chop();
		if(H_(i,nb(2,0))>0){ linestyle = "solid"; }
		else { linestyle = "dashed"; }
		/*y-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth=1pt,linecolor="+color+",linestyle="+linestyle);
	}

	for(unsigned int i(0);i<n_;i++) {
		xy0 = get_pos_in_lattice(i);
		set_pos_LxLy(xy0);
		//set_in_basis(xy0);
		xy0 = (LxLy_*xy0).chop();
		xy0 = inv_e*xy0;
		ps.put(xy0(0)-0.20,xy0(1)+0.15,my::tostring(i));
	}

	Vector<double> Lx(2);
	Lx(0) = LxLy_(0,0);
	Lx(1) = LxLy_(1,0);
	Lx = inv_e*Lx;
	Vector<double> Ly(2);
	Ly(0) = LxLy_(0,1);
	Ly(1) = LxLy_(1,1);
	Ly = inv_e*Ly;

	Matrix<double> polygon(4,2);
	polygon(0,0)=0;
	polygon(0,1)=0;
	polygon(1,0)=Lx(0);
	polygon(1,1)=Lx(1);
	polygon(2,0)=Lx(0)+Ly(0);
	polygon(2,1)=Lx(1)+Ly(1);
	polygon(3,0)=Ly(0);
	polygon(3,1)=Ly(1);
	for(unsigned int i(0);i<polygon.row();i++){ polygon(i,0) -= 1; }
	ps.polygon(polygon,"linecolor=green");

	Vector<double> a(2);
	a(0) = ab_(0,0);
	a(1) = ab_(1,0);
	a = inv_e*a;
	Vector<double> b(2);
	b(0) = ab_(0,1);
	b(1) = ab_(1,1);
	b = inv_e*b;

	polygon(0,0)=0;
	polygon(0,1)=0;
	polygon(1,0)=a(0);
	polygon(1,1)=a(1);
	polygon(2,0)=a(0)+b(0);
	polygon(2,1)=a(1)+b(1);
	polygon(3,0)=b(0);
	polygon(3,1)=b(1);
	for(unsigned int i(0);i<polygon.row();i++){ polygon(i,0) -= 1; }
	ps.polygon(polygon,"linecolor=blue");

	ps.end(true,true,true);
}
/*}*/
