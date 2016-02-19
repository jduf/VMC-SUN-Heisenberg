#include "SquareJastrow.hpp"

SquareJastrow::SquareJastrow(System const& s, Matrix<double> const& nu):
	System(s),
	Square<double>(set_ab(),2,"square-jastrow")
{
	if(status_==3){ init_lattice(); }
	init_bosonic(z_,nu);
	compute_nn();
	compute_sublattice();
	compute_omega_cc();

	system_info_.text("Staggered magnetic field, Becca's idea to mimic an on site chemical potential");
	//std::cout<<"check everything "<<status_<<std::endl;
}

/*{method needed for running*/
void SquareJastrow::create(){ status_--; }

void SquareJastrow::compute_nn(){
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		nn_(obs_[0](i,0),obs_[0](i,1)) = obs_[0](i,1);
	}
	//unsigned int l(z_);
	//for(unsigned int j(0);j<z_;j++){
	//nb = get_neighbourg(nn_(i,j));
	//for(unsigned int k(j);k<j+2;k++){
	//nn_(i,l) = nb(k%z_,0);
	//l++;
	//}
	//}
}

void SquareJastrow::compute_sublattice(){
	for(unsigned int i(0);i<n_;i++){
		if(obs_[0](i,5)<2){ sl_(i) = obs_[0](i,5); }
		else{ std::cout<<"void SquareJastrow::compute_sublattice() : unknown sublatttice."; }
	}
	//std::cout<<"sublattice:"<<std::endl;
	//std::cout<<sl_<<std::endl;
}

void SquareJastrow::compute_omega_cc(){
	if(N_==2){
		omega_(1,1) = -1.0;
		cc_(0,0) = 0;//|up,up>
		cc_(0,1) = 1;//|up,down>
		cc_(1,0) = 1;//|down,up>
		cc_(1,1) = 1;//|down,down>
	}
	/*!\warning omega might need to be complex*/
	//if(N_==3){
	//omega_(1,1) = std::polar(1.0,2.0*M_PI/3.0);
	//omega_(2,2) = std::polar(1.0,2.0*M_PI/3.0);
	//omega_(1,2) = std::polar(1.0,4.0*M_PI/3.0);
	//omega_(2,1) = std::polar(1.0,4.0*M_PI/3.0);
	//cc_(0,0) = 0;
	//cc_(0,1) = 1;
	//cc_(0,2) = 2;
	//cc_(1,0) = 1;
	//cc_(1,1) = 3;
	//cc_(1,2) = 4;
	//cc_(2,0) = 2;
	//cc_(2,1) = 4;
	//cc_(2,2) = 4;
	//}
}

void SquareJastrow::save_param(IOFiles& w) const {
	w.write("nn (nearst neighbours)",nn_);
	w.write("cc (to match nu and x)",cc_);
	w.write("sl (sublattice)",sl_);
	w.write("omega (omega)",omega_);
	GenericSystem<double>::save_param(w);
}

Matrix<double> SquareJastrow::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 2;
	tmp(1,0) = 0;
	tmp(0,1) = 0;
	tmp(1,1) = 1;
	return tmp;
}

unsigned int SquareJastrow::unit_cell_index(Vector<double> const& x) const {
	Vector<double> match(2,0);
	if(my::are_equal(x,match)){ return 0; }
	match(0) = 0.5;
	match(1) = 0;
	if(my::are_equal(x,match)){ return 1; }
	return 2;
}
/*}*/

/*{method needed for checking*/
void SquareJastrow::display_results(){
	//PSTricks ps("./",filename_+"-lattice");
	//ps.add("\\begin{pspicture}(15,15)%"+filename_+"-lattice");
	//Matrix<int> nb;
	//double x0, y0, x1, y1;
	//std::string color;
	//for(unsigned int i(0);i<n_;i++){
	//nb = get_neighbourg(i);
	//x0 = i%Lx_;
	//y0 = i/Ly_;
	//y1 = nb(0,0)/Ly_;
	//if((i+1) % Lx_ ){
	//x1 = nb(0,0)%Lx_;
	//color = "black";
	//} else {
	//x1 = x0 + 1;
	//color = "blue";
	//}
	//ps.line("-", x0, y0, x1, y1 , "linewidth=1pt,linecolor="+color);
	//
	//x1 = nb(1,0)%Lx_;
	//if( i+Lx_<this->n_){
	//y1 = nb(1,0)/Ly_;
	//color = "black";
	//} else {
	//y1 = y0 + 1;
	//color = "blue";
	//}
	//ps.line("-", x0, y0, x1, y1, "linewidth=1pt,linecolor="+color);
	//}
	//
	//double r(0.2);
	//Vector<double> pie(N_);
	//double m(lat.max());
	//for(unsigned int i(0);i<n_;i++){
	//for(unsigned int j(0);j<N_;j++){
	//pie(j) = lat(i,j)/m;
	//if(pie(j) < 1e-4){ pie(j) = 0; }
	//}
	//ps.add("\\rput("+my::tostring(i%Lx_)+","+my::tostring(i/Ly_)+"){%");
	//ps.pie(pie,r,"chartColor=color,userColor={red,blue}");
	//ps.add("}");
	//
	//switch(sl_(i)){
	//case 0: { color = "red"; } break;
	//case 1: { color = "blue"; } break;
	//case 2: { color = "green"; } break;
	//}
	//ps.put(i%Lx_+r*1.2, i/Ly_+r*1.2, "\\tiny{"+my::tostring(i)+"}");
	//}
	//
	//ps.frame(-0.5,-0.5,Lx_-0.5,Ly_-0.5,"linecolor=red");
	//ps.add("\\end{pspicture}");
}

void SquareJastrow::check(){
	std::cout<<"void SquareJastrow::check() : nothing to do"<<std::endl;
}
/*}*/
