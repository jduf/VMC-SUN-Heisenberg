#include "SquarePiFlux.hpp"

SquarePiFlux::SquarePiFlux(System const& s):
	System(s),
	Square<std::complex<double> >(set_ab(),2,"square-piflux")
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		system_info_.text("SquarePiFlux :");
		system_info_.item("Each color has the same Hamiltonian.");
		system_info_.item("There is a flux of "+RST::math("\\pi") + "per square plaquette.");
	}
}

/*{method needed for running*/
void SquarePiFlux::compute_H(){
	double phi(M_PI/4.0);
	H_.set(n_,n_,0);
	unsigned int s0;
	unsigned int s1;
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		if(obs_[0](i,3)){ H_(s0,s1) = std::polar(double(obs_[0](i,4)?bc_:1),obs_[0](i,5)?-phi:phi); }
		else{ H_(s0,s1) = (obs_[0](i,4)?bc_:1); }
	}
	H_ += H_.conjugate_transpose();
}

void SquarePiFlux::create(){
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

Matrix<double> SquarePiFlux::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 2;
	tmp(1,0) = 0;
	tmp(0,1) = 0;
	tmp(1,1) = 1;
	return tmp;
}

unsigned int SquarePiFlux::match_pos_in_ab(Vector<double> const& x) const {
	Vector<double> match(2,0);
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 0; }
	match(0) = 0.5;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 1; }
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<" | "<<x.size()<<std::endl;
	return 2;
}
/*}*/

/*{method needed for checking*/
void SquarePiFlux::display_results(){
	compute_H();

	std::string color("black");
	std::string linestyle("solid");
	std::string arrow("-");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	std::complex<double> t;
	PSTricks ps(info_+path_+dir_,filename_);
	ps.begin(-20,-20,20,20,filename_);
	ps.polygon(cluster_vertex_,"linecolor=green");
	ps.polygon(draw_unit_cell(),"linecolor=black");
	ps.linked_lines("-",draw_boundary(false),"linecolor=yellow");

	double phi(-M_PI/4.0);
	unsigned int s0;
	unsigned int s1;
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		xy0 = x_[s0];

		s1 = obs_[0](i,1);
		xy1 = x_[s1];

		t = H_(s0,s1);
		if(std::abs(t)>1e-4){
			if((xy0-xy1).norm_squared()>1.0001){
				linestyle = "dashed";
				xy1 = (xy0+dir_nn_[obs_[0](i,3)]).chop();
				ps.put(xy1(0)-0.20,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
			} else { linestyle = "solid"; }

			if(t.real()>0){ color = "blue"; }
			else { color = "red"; }

			arrow = "-";
			if(std::arg(t)>0){ arrow = "->"; }
			if(std::arg(t)<0){ arrow = "<-"; }

			ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1), "linewidth=1mm,linecolor="+color+",linestyle="+linestyle);
		}
		if(i%2){
			ps.put(xy0(0)-0.20,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}");
			ps.put(my::chop((2*xy0(0)+1)/2.0),my::chop((xy0(1)+xy1(1))/2.0),"\\tiny{"+my::tostring(2.*(phi-std::arg(t))/M_PI)+"}");
			phi =  std::arg(t);
		}
	}
	ps.end(true,true,true);
}

void SquarePiFlux::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="square-piflux";
	display_results();
}
/*}*/

/*{method needed for analysing*/
std::string SquarePiFlux::extract_level_7(){
	rst_file_ = new RSTFile(info_+path_+dir_,filename_);

	unsigned int nruns;
	unsigned int tmax;

	(*read_)>>nruns>>tmax;
	(*data_write_)<<"% E dE 0|1"<<IOFiles::endl;
	obs_.push_back(Observable(*read_));
	obs_.push_back(Observable(*read_));
	(*data_write_)<<obs_[0][0]<<IOFiles::endl;
	jd_write_->write("E",obs_[0][0]);

	rst_file_->text(read_->get_header());
	rst_file_->save(false,true);
	delete rst_file_;
	rst_file_ = NULL;

	return filename_;
}

std::string SquarePiFlux::extract_level_3(){
	(*read_)>>obs_[0][0];
	(*data_write_)<<n_<<" "<<obs_[0][0]<<" "<<bc_<<IOFiles::endl;

	return filename_;
}
/*}*/
