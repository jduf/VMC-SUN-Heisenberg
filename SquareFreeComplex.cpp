#include "SquareFreeComplex.hpp"

SquareFreeComplex::SquareFreeComplex(System const& s, Vector<double> const& t, Vector<double> const& mu, Vector<double> const& phi):
	System(s,3),
	Square<std::complex<double> >(set_ab(),(N_/m_==2?2:0),"square-free-flux"),
	t_(t),
	mu_(mu),
	phi_(phi)
{
	if(status_==2){
		init_fermionic();

		system_info_.text("FreeComplex : all colors experience the same Hamiltonian");
	}
}

/*{method needed for running*/
void SquareFreeComplex::compute_H(unsigned int const& c){
	H_.set(n_,n_,0);
	Matrix<int> nb;
	unsigned int s(0);
	for(unsigned int i(0);i<n_;i++){
		s = get_site_in_ab(i);
		if(c%2){ if( s==0 || s==3 ){ H_(i,i) = mu_(0)/2.0; } }
		else { if( s==2 || s==1 ){ H_(i,i) = mu_(0)/2.0; } }
		nb = get_neighbourg(i);
		switch(s){ 
			case 0:
				{
					H_(i,nb(0,0)) = std::polar(t_(0)*nb(0,1),phi_(0));
					H_(i,nb(1,0)) = std::polar(t_(1)*nb(1,1),-phi_(0));
				}break;
			case 1:
				{
					H_(i,nb(0,0)) = std::polar(1.0*nb(0,1),-phi_(0));
					H_(i,nb(1,0)) = std::polar(t_(1)*nb(1,1),phi_(0));
				}break;
			case 2:
				{
					H_(i,nb(0,0)) = std::polar(t_(0)*nb(0,1),-phi_(0));
					H_(i,nb(1,0)) = std::polar(1.0*nb(1,1),phi_(0));
				}break;
			case 3:
				{
					H_(i,nb(0,0)) = std::polar(1.0*nb(0,1),phi_(0));
					H_(i,nb(1,0)) = std::polar(1.0*nb(1,1),-phi_(0));
				}break;
			default:{ std::cerr<<"void SquareFreeComplex::compute_H(unsigned int const& c) : undefined site in unit cell"<<std::endl; }break;
		}
	}
	H_ += H_.trans_conj(); 
}

void SquareFreeComplex::create(){
	for(unsigned int c(0);c<N_;c++){
		compute_H(c);
		diagonalize(true);
		if(status_==1){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = H_(i,j);
				}
			}
			if(c!=N_-1){ status_++;}
		}
	}
}

unsigned int SquareFreeComplex::match_pos_in_ab(Vector<double> const& x) const{
	Vector<double> match(2,0);
	if(my::are_equal(x,match)){ return 0; }
	match(0) = 0.5;
	match(1) = 0;
	if(my::are_equal(x,match)){ return 1; }
	match(0) = 0;
	match(1) = 0.5;
	if(my::are_equal(x,match)){ return 2; }
	match(0) = 0.5;
	match(1) = 0.5;
	if(my::are_equal(x,match)){ return 3; }
	return 4;
}

Matrix<double> SquareFreeComplex::set_ab(){
	Matrix<double> tmp(2,2);
	tmp(0,0) = 2;
	tmp(1,0) = 0;
	tmp(0,1) = 0;
	tmp(1,1) = 2;
	return tmp;
}
/*}*/

/*{method needed for checking*/
void SquareFreeComplex::lattice(std::string const& path){
	compute_H(1);
	Matrix<int> nb;
	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	std::string arrow("-");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(path,"lattice");
	ps.add("\\begin{pspicture}(-20,-20)(20,20)%"+filename_);
	for(unsigned int i(0);i<n_;i++) {
		xy0 = get_pos_in_lattice(i);
		set_pos_LxLy(xy0);
		set_in_basis(xy0);
		xy0 = (LxLy_*xy0).chop();
		ps.put(xy0(0)-0.20,xy0(1)+0.15,my::tostring(i));
		nb = get_neighbourg(i);

		if(nb(0,1)<0){ color = "red"; }
		else { color = "black"; }
		xy1 = get_pos_in_lattice(nb(0,0));
		set_pos_LxLy(xy1);
		set_in_basis(xy1);
		xy1 = LxLy_*xy1;
		xy1 = xy1.chop();
		if( xy1(0)<xy0(0) || i+1==xloop_){ 
			xy1 = xy0;
			xy1(0) += 1.0;
			ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(nb(0,0)));
			linestyle="dashed";
		} else { linestyle="solid"; }

		if(arg(H_(i,nb(0,0)))>0){ arrow = "->"; }
		else { arrow = "<-"; }
		linewidth = my::tostring(std::abs(H_(i,nb(0,0))))+"pt";
		/*x-link*/ ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);

		if(nb(1,1)<0){ color = "red"; }
		else { color = "black"; }
		xy1 = get_pos_in_lattice(nb(1,0));
		set_pos_LxLy(xy1);
		set_in_basis(xy1);
		xy1 = (LxLy_*xy1).chop();
		if( xy1(1)<xy0(1) ){ 
			xy1 = xy0;
			xy1(1) += 1.0;
			ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(nb(1,0)));
			linestyle="dashed";
		} else { linestyle="solid"; }
		if(arg(H_(i,nb(1,0)))>0){ arrow = "->"; }
		else { arrow = "<-"; }
		linewidth = my::tostring(std::abs(H_(i,nb(1,0))))+"pt";
		/*y-link*/ ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);

		if(my::real(H_(i,i))){ ps.circle(xy0,my::real(H_(i,i)),"linecolor=magenta,fillstyle=solid,fillcolor=magenta"); }
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

void SquareFreeComplex::check(){
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
