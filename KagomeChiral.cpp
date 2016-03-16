#include "KagomeChiral.hpp"

KagomeChiral::KagomeChiral(System const& s, double const& phi):
	System(s),
	Kagome<std::complex<double> >(set_ab(),6,"kagome-chiral"),
	phi_(phi)
{
	if(phi<=3.0){
		if(status_==3){ init_lattice(); }
		if(status_==2){
			init_fermionic();

			system_info_.text("KagomeChiral :");
			system_info_.item("Each color has the same Hamiltonian.");
			system_info_.item("6 sites per unit cell.");
			system_info_.item("Flux of "+RST::math(my::tostring(phi)+"\\pi/3")+" per plaquette.");
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : the flux per plaquette shouldn't be bigger than pi (phi<=3)"<<std::endl; }
}

/*{method needed for running*/
void KagomeChiral::compute_H(){
	H_.set(n_,n_,0);

	double t(-1.0);
	double phi(phi_*M_PI/3.0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		switch(obs_[0](i,5)){
			case 0:{ H_(obs_[0](i,0),obs_[0](i,1)) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)==0?-phi:0.0)); }break;
			case 1:{ H_(obs_[0](i,0),obs_[0](i,1)) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)==2?phi:0.0)); }break;
			case 2:{ H_(obs_[0](i,0),obs_[0](i,1)) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)==0?3.0*phi:0.0)); }break;
			case 3:{ H_(obs_[0](i,0),obs_[0](i,1)) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)==0?-phi:0.0)); }break;
			case 4:{ H_(obs_[0](i,0),obs_[0](i,1)) = std::polar((obs_[0](i,4)?bc_*t:t),(obs_[0](i,3)==2?-2.0*phi:0.0)); }break;
			case 5:{ H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*t:t); }break;
		}
	}
	H_ += H_.conjugate_transpose();
}

void KagomeChiral::create(){
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

void KagomeChiral::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("phi="+my::tostring(phi_));
		Vector<double> param(1,phi_);

		w.add_header()->title(s,'<');
		w<<param;
		w.add_header()->add(system_info_.get());
	} else { w<<phi_<<" "; }
}

Matrix<double> KagomeChiral::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 2.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 0.0;
	tmp(1,1) = 2.0*sqrt(3.0);
	return tmp;
}

unsigned int KagomeChiral::unit_cell_index(Vector<double> const& x) const {
	if(my::are_equal(x(1),0.0,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.0,eq_prec_,eq_prec_)){ return 0; }
		if(my::are_equal(x(0),0.5,eq_prec_,eq_prec_)){ return 2; }
	}
	if(my::are_equal(x(1),0.25,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.25,eq_prec_,eq_prec_)){ return 1; }
	}
	if(my::are_equal(x(1),0.5,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.0,eq_prec_,eq_prec_)){ return 5; }
		if(my::are_equal(x(0),0.5,eq_prec_,eq_prec_)){ return 3; }
	}
	if(my::are_equal(x(1),0.75,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.75,eq_prec_,eq_prec_)){ return 4; }
	}
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void KagomeChiral::lattice(){
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
	Matrix<double> uc(draw_unit_cell(-1.0,-sqrt(3.0)/4.0));
	ps.polygon(uc,"linecolor=black");
	ps.linked_lines("-",draw_boundary(false),"linecolor=yellow");

	unsigned int s0;
	unsigned int s1;
	/*draws only the lattice, shows links and bc*/
	//for(unsigned int i(0);i<obs_[0].nlinks();i++){
		//s0 = obs_[0](i,0);
		//xy0 = x_[s0];
		//s1 = obs_[0](i,1);
		//xy1 = x_[s1];
//
		//if((xy0-xy1).norm_squared()>1.0001){
			//linestyle = "dashed";
			//xy1 = (xy0+dir_nn_[obs_[0](i,3)]*1.2).chop();
			//ps.put(xy1(0),xy1(1),"\\tiny{"+my::tostring(s1)+"}");
			//xy1 = (xy0+dir_nn_[obs_[0](i,3)]).chop();
		//} else { linestyle = "solid"; }
//
		//ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1),"linewidth=1pt,linecolor=black,linestyle="+linestyle);
//
		//if(i%2){
			//switch(obs_[0](i,5)%3){
				//case 0: { ps.put(xy0(0)-0.2,xy0(1)+0.2,"\\tiny{"+my::tostring(s0)+"}"); }break;
				//case 1: { ps.put(xy0(0)-0.2,xy0(1),"\\tiny{"+my::tostring(s0)+"}"); }break;
				//case 2: { ps.put(xy0(0)-0.2,xy0(1)-0.2,"\\tiny{"+my::tostring(s0)+"}"); }break;
			//}
		//}
	//}
	/*draws long range correlations over the lattice*/
	if(o(1)){ draw_long_range_correlation(ps,obs_[o(1)]); }

	Vector<double> shift(2,0.0);
	if(o(0) || o(2)){
		/*unit cell, shows bond energy and color occupation*/
		double be;
		//Vector<double> shift(equivalent_vertex_[0]+equivalent_vertex_[1]);
		//ps.polygon(draw_unit_cell(shift(0)+0.5,shift(1)+0.5),"linecolor=black");
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
	}

	/*unit cell, shows hopping amplitude, chemical potential and fluxes*/
	std::complex<double> t;
	double flux;
	double sign;
	unsigned long long a;
	unsigned long long b;
	std::string arrow("-");
	//shift = equivalent_vertex_[0]+equivalent_vertex_[2];
	//ps.polygon(draw_unit_cell(shift(0)+0.5,shift(1)+0.5),"linecolor=black");
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		xy0 = x_[s0];
		s1 = obs_[0](i,1);
		xy1 = x_[s1];

		if((xy0-xy1).norm_squared()>1.0001){
			linestyle = "dashed";
			xy1 = (xy0+dir_nn_[obs_[0](i,3)]*1.2).chop();
			ps.put(xy1(0),xy1(1),"\\tiny{"+my::tostring(s1)+"}");
			xy1 = (xy0+dir_nn_[obs_[0](i,3)]).chop();
		} else { linestyle = "solid"; }

		//if(my::in_polygon(uc.row(),uc.ptr(),uc.ptr()+uc.row(),xy0(0),xy0(1)))
		{
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
				ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
			}

			if(i%2){
				switch(obs_[0](i,5)%3){
					case 0: { ps.put(xy0(0)-0.2,xy0(1)+0.2,"\\tiny{"+my::tostring(s0)+"}"); }break;
					case 1: { ps.put(xy0(0)-0.2,xy0(1),"\\tiny{"+my::tostring(s0)+"}"); }break;
					case 2: { ps.put(xy0(0)-0.2,xy0(1)-0.2,"\\tiny{"+my::tostring(s0)+"}"); }break;
				}
			}

			switch(obs_[0](i,5)%3){
				case 0:
					{
						if(obs_[0](i,3)==0){
							flux = 0.0;
							unsigned int j(0);
							do {
								xy0 += dir_nn_[2*j];
								s1 = site_index(xy0);
								flux += std::arg(-H_(s0,s1));
								s0 = s1;
							} while (++j<3);
							flux = my::chop(flux/M_PI);
							if(!my::are_equal(flux,0.0,eq_prec_,eq_prec_)){
								if(my::to_fraction(flux,a,b,sign) && b!=1){
									ps.put((xy0(0)+xy1(0))/2.0,xy0(1)+sqrt(3.0)/4.0,"\\tiny{"+std::string(sign<0?"-":"")+"$\\frac{"+my::tostring(a)+"}{"+my::tostring(b)+"}$}");
								} else {
									ps.put((xy0(0)+xy1(0))/2.0,xy0(1)+sqrt(3.0)/4.0,"\\tiny{"+my::tostring(flux)+"}");
								}
							}
						}
					}break;
				case 1:
					{
						if(obs_[0](i,3)==1){
							flux = 0.0;
							unsigned int j(0);
							do {
								xy0 += dir_nn_[2*j+1];
								s1 = site_index(xy0);
								std::cout<<s0<<" "<<s1<<std::endl;
								flux += std::arg(-H_(s0,s1));
								s0 = s1;
							} while (++j<3);
							std::cout<<std::endl;
							flux = my::chop(flux/M_PI);
							if(!my::are_equal(flux,0.0,eq_prec_,eq_prec_)){
								if(my::to_fraction(flux,a,b,sign) && b!=1){
									ps.put(xy0(0),(xy0(1)+xy1(1))/2.0,"\\tiny{"+std::string(sign<0?"-":"")+"$\\frac{"+my::tostring(a)+"}{"+my::tostring(b)+"}$}");
								} else {
									ps.put(xy0(0),(xy0(1)+xy1(1))/2.0,"\\tiny{"+my::tostring(flux)+"}");
								}
							}
						}
					}break;
				case 2:
					{
						if(obs_[0](i,3)==0){
							flux = 0.0;
							unsigned int j(0);
							do {
								xy0 += dir_nn_[j];
								s1 = site_index(xy0);
								flux += std::arg(-H_(s0,s1));
								s0 = s1;
							} while (++j<6);
							flux = my::chop(flux/M_PI);
							if(!my::are_equal(flux,0.0,eq_prec_,eq_prec_)){
								if(my::to_fraction(flux,a,b,sign) && b!=1){
									ps.put((xy0(0)+xy1(0))/2.0,xy0(1)+1.0,"\\tiny{"+std::string(sign<0?"-":"")+"$\\frac{"+my::tostring(a)+"}{"+my::tostring(b)+"}$}");
								} else {
									ps.put((xy0(0)+xy1(0))/2.0,xy0(1)+1.0,"\\tiny{"+my::tostring(flux)+"}");
								}
							}
						}
					}break;
			}
		}
	}
	ps.end(true,true,true);
}

void KagomeChiral::display_results(){
	lattice();

	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title(RST::math("\\phi")+"="+my::tostring(phi_));
		std::string run_cmd("./mc -s:wf kagome-chiral");
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

void KagomeChiral::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="kagome-chiral";
	//display_results();

	compute_H();
	plot_band_structure();
}
/*}*/
