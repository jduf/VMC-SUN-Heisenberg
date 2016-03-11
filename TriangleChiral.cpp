#include "TriangleChiral.hpp"

TriangleChiral::TriangleChiral(System const& s, double const& phi):
	System(s),
	Triangle<std::complex<double> >(set_ab(),3,"triangle-chiral"),
	phi_(phi)
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		system_info_.text("TriangleChiral :");
		system_info_.item("Each color has the same Hamiltonian.");
		system_info_.item("There is a flux of "+RST::math(my::tostring(phi)+"\\pi/3") + "per triangular plaquette.");

		filename_ += "-phi"+my::tostring(phi_);
	}
}

/*{method needed for running*/
void TriangleChiral::compute_H(){
	H_.set(n_,n_,0);

	double t(-1.0);
	unsigned int s0(0);
	unsigned int s1(0);
	double phi(phi_*M_PI/3.0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);

		switch(obs_[0](i,3)){
			case 0:{ H_(s0,s1) = (obs_[0](i,4)?bc_*t:t); }break;
			case 1:{ H_(s0,s1) = std::polar((obs_[0](i,4)?bc_*t:t),phi*(2.0*obs_[0](i,5)+1)); }break;
			case 2:{ H_(s0,s1) = std::polar((obs_[0](i,4)?bc_*t:t),phi* 2.0*obs_[0](i,5)); }break;
		}
	}
	H_ += H_.conjugate_transpose();
}

void TriangleChiral::create(){
	compute_H();
	diagonalize(true);
	if(status_==1){
		for(unsigned int c(0);c<N_;c++){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = H_(i,j);
				}
			}
		}
	}
}

void TriangleChiral::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("phi="+my::tostring(phi_));
		Vector<double> param(1,phi_);

		w.add_header()->title(s,'<');
		w<<param;
		w.add_header()->add(system_info_.get());
	} else { w<<phi_<<" "; }
}

Matrix<double> TriangleChiral::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 3.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 0.5;
	tmp(1,1) = sqrt(3.0)/2.0;
	return tmp;
}

unsigned int TriangleChiral::unit_cell_index(Vector<double> const& x) const {
	if(my::are_equal(x(0),0.0    ,eq_prec_,eq_prec_)){ return 0; }
	if(my::are_equal(x(0),1.0/3.0,eq_prec_,eq_prec_)){ return 1; }
	if(my::are_equal(x(0),2.0/3.0,eq_prec_,eq_prec_)){ return 2; }
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void TriangleChiral::lattice(){
	Vector<unsigned int> o(3,0);
	for(unsigned int i(1);i<obs_.size();i++){
		switch(obs_[i].get_type()){
			case 1:{ o(0)=i; }break;//bond energy
			case 2:{ o(1)=i; }break;//long range correlation
			case 3:{ o(2)=i; }break;//color occupation
		}
	}
	compute_H();

	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	std::string arrow("-");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(info_+path_+dir_,filename_);
	ps.begin(-20,-20,20,20,filename_);
	ps.polygon(cluster_vertex_,"linecolor=green");
	Matrix<double> uc(draw_unit_cell());
	ps.polygon(uc,"linecolor=black");
	ps.linked_lines("-",draw_boundary(false),"linecolor=yellow");

	std::complex<double> t;
	unsigned int s0;
	unsigned int s1;
	Matrix<int> nb;
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		xy0 = x_[s0];

		s1 = obs_[0](i,1);
		xy1 = x_[s1];

		t = H_(s0,s1);
		if(o(0) || o(2)){
			if(my::in_polygon(uc.row(),uc.ptr(),uc.ptr()+uc.row(),xy0(0),xy0(1))){
				if(o(0)){ t = obs_[o(0)][obs_[0](i,2)].get_x(); }
				if(i%2 && o(2)){
					Vector<double> p(N_);
					for(unsigned int j(0);j<N_;j++){ p(j) = obs_[o(2)][j+N_*obs_[0](i,5)].get_x(); }
					ps.pie(xy0(0),xy0(1),p,0.2,"chartColor=color");
				}
			} else if(my::in_polygon(uc.row(),uc.ptr(),uc.ptr()+uc.row(),xy1(0),xy1(1))){ t = 0; }
			linewidth = my::tostring(std::abs(t))+"mm";
		}
		if(std::abs(t)>1e-4){
			if((xy0-xy1).norm_squared()>1.0001){
				linestyle = "dashed";
				xy1 = (xy0+dir_nn_[obs_[0](i,3)]).chop();
				ps.put(xy1(0)+0.2,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
			} else { linestyle = "solid"; }

			if(t.real()>0){ color = "blue"; }
			else          { color = "red"; }

			if(my::are_equal(t.imag(),0)){ arrow = "-"; }
			else { arrow = "->"; }

			ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1),"linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}
		if(i%3==2){ ps.put(xy0(0)+0.2,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}"); }
		if(i%3<2){
			double sign;
			unsigned int j(0);
			unsigned long long a;
			unsigned long long b;
			double flux(0.0);
			do {
				nb = get_neighbourg(s0);
				s1 = nb(i%3+2*j,0);
				flux += std::arg(-H_(s0,s1));
				s0 = s1;
			} while (++j<3);
			flux = my::chop(flux/M_PI);
			if(my::to_fraction(flux,a,b,sign)){
				if(i%3){ ps.put(xy0(0),xy0(1)+sqrt(3.0)/4.0,"\\tiny{"+std::string(sign<0?"-":"")+"$\\frac{"+my::tostring(a)+"}{"+my::tostring(b)+"}$}"); }
				else   { ps.put((xy0(0)+xy1(0))/2.0,xy0(1)+sqrt(3.0)/4.0,"\\tiny{"+std::string(sign<0?"-":"")+"$\\frac{"+my::tostring(a)+"}{"+my::tostring(b)+"}$}"); }
			} else {
				if(i%3){ ps.put(xy0(0),xy0(1)+sqrt(3.0)/4.0,"\\tiny{"+my::tostring(flux)+"}"); }
				else   { ps.put((xy0(0)+xy1(0))/2.0,xy0(1)+sqrt(3.0)/4.0,"\\tiny{"+my::tostring(flux)+"}"); }
			}
		}
	}
	if(o(1)){ draw_long_range_correlation(ps,obs_[o(1)]); }
	ps.end(true,true,true);
}

void TriangleChiral::display_results(){
	lattice();

	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title(RST::math("\\phi")+"="+my::tostring(phi_));
		std::string run_cmd("./mc -s:wf triangle-chiral");
		run_cmd += " -u:N " + my::tostring(N_);
		run_cmd += " -u:m " + my::tostring(m_);
		run_cmd += " -u:n " + my::tostring(n_);
		run_cmd += " -i:bc "+ my::tostring(bc_);
		run_cmd += " -d:phi "+ my::tostring(phi_);;
		run_cmd += " -d -u:tmax 10";

		rst_file_->title(title,'-');
		rst_file_->change_text_onclick("run command",run_cmd);

		rst_file_->figure(dir_+filename_+".png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(dir_+filename_+".pdf")+RST::scale("200"));
	}
}

void TriangleChiral::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="triangle-chiral";
	display_results();

	//compute_H();
	//plot_band_structure();
}
/*}*/
