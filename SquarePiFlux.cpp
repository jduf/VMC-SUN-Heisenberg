#include "SquarePiFlux.hpp"

SquarePiFlux::SquarePiFlux(System const& s):
	System(s),
	Square<std::complex<double> >(set_ab(),2,"square-piflux")
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		/*This was an attempt to recover results obtained by Wang and
		 * Vishwanath. But as they are using another representation (using
		 * Majorana fermions), the variational wavefunction that they provide
		 * (the one implemented in this class) can't be used to obtain good
		 * result on SU(4) m=1. The energy is way to high.*/
		system_info_.text("SquarePiFlux :");
		system_info_.item("Each color has the same Hamiltonian.");
		system_info_.item("There is a flux of "+RST::math("\\pi") + "per square plaquette.");
	}
}

/*{method needed for running*/
void SquarePiFlux::compute_H(){
	H_.set(n_,n_,0);

	double phi(M_PI/2.0);
	unsigned int s0;
	unsigned int s1;
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		if(obs_[0](i,3)){ H_(s0,s1) = std::polar(obs_[0](i,4)?bc_:1.0,phi)/2.0; }
		else{ H_(s0,s1) = std::polar(obs_[0](i,4)?bc_:1.0,obs_[0](i,5)?-phi:phi)/2.0; }
	}
	H_ += H_.conjugate_transpose();
}

void SquarePiFlux::create(){
	compute_H();
	diagonalize(true);
	status_=1;
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

Matrix<double> SquarePiFlux::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 1;
	tmp(1,0) = 0;
	tmp(0,1) = 0;
	tmp(1,1) = 2;
	return tmp;
}

unsigned int SquarePiFlux::unit_cell_index(Vector<double> const& x) const {
	return my::are_equal(x(1),0.5,eq_prec_,eq_prec_);
}
/*}*/

/*{method needed for checking*/
void SquarePiFlux::lattice(){
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
	Matrix<double> uc(draw_unit_cell(0,0.5));
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
		if(obs_.size()>1){
			if(my::in_polygon(uc.row(),uc.ptr(),uc.ptr()+uc.row(),xy0(0),xy0(1))){ 
				t = obs_[1][i%4].get_x(); 
				if(i%2 && obs_.size()>2){
					Vector<double> p(N_);
					for(unsigned int j(0);j<N_;j++){ p(j) = obs_[2][j+N_*obs_[0](i,5)].get_x(); }
					ps.pie(xy0(0),xy0(1),p,0.2,"chartColor=color");
				}
			} else if(my::in_polygon(uc.row(),uc.ptr(),uc.ptr()+uc.row(),xy1(0),xy1(1))){ t = 0; }
		}
		if(std::abs(t)>1e-4){
			if((xy0-xy1).norm_squared()>1.0001){
				linestyle = "dashed";
				xy1 = (xy0+dir_nn_[obs_[0](i,3)]).chop();
				ps.put(xy1(0)+0.2,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
			} else { linestyle = "solid"; }

			if(t.real()>0){ color = "blue"; }
			else          { color = "red"; }

			if(my::are_equal(t.imag(),0.0)){ arrow = "-"; }
			else { arrow = "->"; }
			//ps.put((xy0(0)+xy1(0))/2.0,(xy0(1)+xy1(1))/2.0,"\\tiny{"+my::tostring(my::chop(std::arg(-t)/M_PI))+"}");

			ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}

		if(obs_[0](i,3)){ ps.put(xy0(0)+0.1,(xy0(1)+xy1(1))/2.0,"\\tiny{"+std::string(1,my::int_to_alphabet(obs_[0](i,2),true))+"}"); }
		else            { ps.put((xy0(0)+xy1(0))/2.0,xy0(1)+0.1,"\\tiny{"+std::string(1,my::int_to_alphabet(obs_[0](i,2),true))+"}"); }

		if(i%2){
			ps.put(xy0(0)+0.10,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}"); 
		} else {
			unsigned int j(0);
			double flux(0.0);
			do {
				nb = get_neighbourg(s0);
				s1 = nb(j,0);
				flux += std::arg(-H_(s0,s1));
				s0 = s1;
			} while (++j<4);
			flux = my::chop(flux/M_PI);
			ps.put((xy0(0)+xy1(0))/2.0,xy0(1)+0.5,"\\tiny{"+my::tostring(flux)+"}");
		}
	}
	ps.end(true,true,true);
}

void SquarePiFlux::display_results(){
	lattice();

	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title("pi-flux");
		std::string run_cmd("./mc -s:wf square-piflux");
		run_cmd += " -u:N " + my::tostring(N_);
		run_cmd += " -u:m " + my::tostring(m_);
		run_cmd += " -u:n " + my::tostring(n_);
		run_cmd += " -i:bc "+ my::tostring(bc_);
		run_cmd += " -d -u:tmax 10";

		rst_file_->title(title,'-');
		rst_file_->change_text_onclick("run command",run_cmd);

		rst_file_->figure(dir_+filename_+".png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(dir_+filename_+".pdf")+RST::scale("200"));
	}
}

void SquarePiFlux::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="square-piflux";
	display_results();

	plot_band_structure();
}
/*}*/
