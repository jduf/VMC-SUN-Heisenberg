#include "SquareChiral.hpp"

SquareChiral::SquareChiral(System const& s, double const& phi):
	System(s),
	Square<std::complex<double> >(set_ab(ref_(3),N_/m_),N_/m_,"square-chiral"),
	phi_(phi)
{
	if(2*phi<=N_/m_){
		if(status_==3){ init_lattice(); }
		if(status_==2){
			init_fermionic();

			system_info_.text("SquareChiral :");
			system_info_.item("Each color has the same Hamiltonian.");
			system_info_.item("Flux of "+RST::math(my::tostring(phi)+"\\times 2\\pi/"+my::tostring(N_/m_))+ " per plaquette.");

			filename_ += "-phi"+my::tostring(phi_);
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : the flux per plaquette shouldn't be bigger than pi"<<std::endl; }
}

/*{method needed for running*/
void SquareChiral::compute_H(){
	H_.set(n_,n_,0);

	double t(-1.0);
	unsigned int s0(0);
	unsigned int s1(0);
	double phi(phi_*2.0*m_*M_PI/N_);
	switch(spuc_){
		case 4:
			{
				for(unsigned int i(0);i<obs_[0].nlinks();i++){
					s0 = obs_[0](i,0);
					s1 = obs_[0](i,1);
					switch(obs_[0](i,5)){
						case 0: { H_(s0,s1) = (obs_[0](i,4)?bc_*t:t); }break;
						case 1: { H_(s0,s1) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)?2.0*phi:0.0)); }break;
						case 2: { H_(s0,s1) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)?0.0:phi)); }break;
						case 3: { H_(s0,s1) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)?0.0:phi)); }break;
					}
				}
			}break;
		case 5:
			{
				for(unsigned int i(0);i<obs_[0].nlinks();i++){
					s0 = obs_[0](i,0);
					s1 = obs_[0](i,1);
					switch(obs_[0](i,5)){
						case 0: { H_(s0,s1) = (obs_[0](i,4)?bc_*t:t); }break;
						case 1: { H_(s0,s1) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)?-phi:phi)); }break;
						case 2: { H_(s0,s1) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)?phi:0.0)); }break;
						case 3: { H_(s0,s1) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)?0.0:-phi)); }break;
						case 4: { H_(s0,s1) = (obs_[0](i,4)?bc_*t:t); }break;
					}
				}
			}break;
	}
	H_ += H_.conjugate_transpose();
}

void SquareChiral::create(){
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

void SquareChiral::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("phi="+my::tostring(phi_));
		Vector<double> param(1,phi_);

		w.add_header()->title(s,'<');
		w<<param;
		w.add_header()->add(system_info_.get());
	} else { w<<phi_<<" "; }
}

Matrix<double> SquareChiral::set_ab(unsigned int const& ref3, unsigned int const& k) const {
	Matrix<double> tmp;
	switch(k){
		case 2:
			{
				tmp.set(2,2);
				tmp(0,0) = 1.0;
				tmp(1,0) =-1.0;
				tmp(0,1) = 1.0;
				tmp(1,1) = 1.0;
			}break;
		case 3:
			{
				tmp.set(2,2);
				tmp(0,0) = 3.0;
				tmp(1,0) = 0.0;
				tmp(0,1) = 0.0;
				tmp(1,1) = 1.0;
			}break;
		case 4:
			{
				tmp.set(2,2);
				tmp(0,0) = 2.0;
				tmp(1,0) = 0.0;
				tmp(0,1) = 0.0;
				tmp(1,1) = 2.0;
			}break;
		case 5:
			{
				tmp.set(2,2);
				if(ref3==2){
					tmp(0,0) = 2.0;
					tmp(1,0) = 1.0;
					tmp(0,1) =-1.0;
					tmp(1,1) = 2.0;
				} else {
					tmp(0,0) = 2.0;
					tmp(1,0) =-1.0;
					tmp(0,1) = 1.0;
					tmp(1,1) = 2.0;
				}
			}break;
		default:
			{ std::cerr<<__PRETTY_FUNCTION__<<" : unknown unit cell"<<std::endl; }
	}
	return tmp;
}

unsigned int SquareChiral::unit_cell_index(Vector<double> const& x) const {
	Vector<double> match(2,0);
	switch(N_/m_){
		case 4:
			{
				unsigned int i(0);
				if(my::are_equal(x(0),0.5,eq_prec_,eq_prec_)){ i+=1; }
				if(my::are_equal(x(1),0.5,eq_prec_,eq_prec_)){ i+=2; }
				return i;
			}
		case 5:
			{
				if(ref_(3)==2){
					if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 0; }
					match(0) = 0.2;
					match(1) = 0.4;
					if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 1; }
					match(0) = 0.6;
					match(1) = 0.2;
					if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 2; }
					match(0) = 0.4;
					match(1) = 0.8;
					if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 3; }
					match(0) = 0.8;
					match(1) = 0.6;
					if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 4; }
				} else {
					if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 0; }
					match(0) = 0.4;
					match(1) = 0.2;
					if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 1; }
					match(0) = 0.8;
					match(1) = 0.4;
					if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 2; }
					match(0) = 0.2;
					match(1) = 0.6;
					if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 3; }
					match(0) = 0.6;
					match(1) = 0.8;
					if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 4; }
				}
			}
	}
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return 5;
}
/*}*/

/*{method needed for checking*/
void SquareChiral::lattice(){
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
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(info_+path_+dir_,filename_);
	ps.begin(-20,-20,20,20,filename_);
	ps.polygon(cluster_vertex_,"linecolor=green");
	Matrix<double> uc(draw_unit_cell(0.5,0.5));
	ps.polygon(uc,"linecolor=black");
	ps.linked_lines("-",draw_boundary(false),"linecolor=yellow");

	unsigned int s0;
	unsigned int s1;
	/*draws only the lattice, shows links and bc*/
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		xy0 = x_[s0];

		s1 = obs_[0](i,1);
		xy1 = x_[s1];

		if((xy0-xy1).norm_squared()>1.0001){
			linestyle = "dashed";
			xy1 = (xy0+dir_nn_[obs_[0](i,3)]).chop();
			ps.put(xy1(0)+0.2,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
		} else { linestyle = "solid"; }

		ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1),"linewidth=1pt,linecolor=black,linestyle="+linestyle);

		if(i%2){ ps.put(xy0(0)+0.2,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}"); }
	}
	/*draws long range correlations over the lattice*/
	if(o(1)){ draw_long_range_correlation(ps,obs_[o(1)]); }

	/*unit cell, shows bond energy and color occupation*/
	double be;
	Vector<double> shift(equivalent_vertex_[0]+equivalent_vertex_[1]);
	ps.polygon(draw_unit_cell(shift(0)+0.5,shift(1)+0.5),"linecolor=black");
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		xy0 = x_[obs_[0](i,0)];
		xy1 = x_[obs_[0](i,1)];

		if(my::in_polygon(uc.row(),uc.ptr(),uc.ptr()+uc.row(),xy0(0),xy0(1))){
			xy0 += shift;
			xy1 += shift;
			if(o(0)){
				be = obs_[o(0)][obs_[0](i,2)].get_x();
				linewidth = my::tostring(std::abs(be))+"mm";
				if(std::abs(be)>1e-4){
					if(be>0){ color = "blue"; }
					else    { color = "red"; }
					ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle=solid");
				}
			}
			if(i%2 && o(2)){
				Vector<double> p(N_);
				for(unsigned int j(0);j<N_;j++){ p(j) = obs_[o(2)][j+N_*obs_[0](i,5)].get_x(); }
				ps.pie(xy0(0),xy0(1),p,0.2,"chartColor=color");
			}
		}
	}

	/*unit cell, shows hopping amplitude, chemical potential and fluxes*/
	std::complex<double> t;
	double flux;
	double sign;
	unsigned long long a;
	unsigned long long b;
	std::string arrow("-");
	shift = equivalent_vertex_[0]+equivalent_vertex_[2];
	ps.polygon(draw_unit_cell(shift(0)+0.5,shift(1)+0.5),"linecolor=black");
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		xy0 = x_[s0];

		s1 = obs_[0](i,1);
		xy1 = x_[s1];

		if(my::in_polygon(uc.row(),uc.ptr(),uc.ptr()+uc.row(),xy0(0),xy0(1))){
			xy0 += shift;
			xy1 += shift;
			t = H_(s0,s1);
			if(std::abs(t)>1e-4){
				linewidth = my::tostring(std::abs(t))+"mm";

				if(my::real(t)>0){ color = "blue"; }
				else             { color = "red"; }

				if(my::are_equal(my::imag(t),0)){ arrow = "-"; }
				else { 
					if(my::imag(t)>0){ arrow = "->"; }
					else { arrow = "<-"; }
				}
				ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle=solid");
			}

			if(!(i%2)){
				unsigned int j(0);
				xy0 = x_[s0];
				flux = 0;
				do {
					xy0 += dir_nn_[j];
					s1 = site_index(xy0);
					flux += std::arg(-H_(s0,s1));
					s0 = s1;
				} while (++j<4);
				flux = my::chop(flux/M_PI);
				if(!my::are_equal(flux,0.0,eq_prec_,eq_prec_)){
					if(my::to_fraction(flux,a,b,sign) && b!=1){
						ps.put(xy0(0),xy0(1)+sqrt(3.0)/4.0,"\\tiny{"+std::string(sign<0?"-":"")+"$\\frac{"+my::tostring(a)+"}{"+my::tostring(b)+"}$}");
					} else {
						ps.put(xy0(0),xy0(1)+sqrt(3.0)/4.0,"\\tiny{"+my::tostring(flux)+"}");
					}
				}
			}
		}
	}
	ps.end(true,true,true);
}

void SquareChiral::display_results(){
	lattice();

	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title(RST::math("\\phi")+"="+my::tostring(phi_));
		std::string run_cmd("./mc -s:wf square-chiral");
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

void SquareChiral::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="square-chiral";
	display_results();

	//compute_H();
	//plot_band_structure();
}
/*}*/
