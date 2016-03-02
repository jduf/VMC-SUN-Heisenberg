#include "SquareMu.hpp"

SquareMu::SquareMu(System const& s, double const& mu):
	System(s),
	Square<double>(set_ab(ref_(3)),5,"square-mu"),
	mu_(mu)
{
	if(status_==3){ init_lattice(); }
	if(status_==2){ 
		init_fermionic();
		same_wf_ = false;

		system_info_.text("SquareMu :");
		system_info_.item("Each color has a different Hamiltonian.");

		filename_ += "-mu"+std::string(mu_>=0?"+":"")+my::tostring(mu_);
	}
}

/*{method needed for running*/
void SquareMu::compute_H(unsigned int const& c){
	H_.set(n_,n_,0);

	double t(-1.0);
	unsigned int s0(0);
	unsigned int s1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		H_(s0,s1) = (obs_[0](i,4)?bc_*t:t);
		if((unsigned int)(obs_[0](i,5))==c%spuc_){ H_(s0,s0) = mu_/2; }
	}
	H_ += H_.transpose();
}

void SquareMu::create(){
	for(unsigned int c(0);c<N_;c++){
		status_ = 2;
		compute_H(c);
		diagonalize(true);
		if(status_==1){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = H_(i,j);
				}
			}
		}
	}
}

void SquareMu::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("mu=("+my::tostring(mu_)+")");
		Vector<double> param(1,mu_);

		w.add_header()->title(s,'<');
		w<<param;
		GenericSystem<double>::save_param(w);
	} else { w<<mu_<<" "; }
}

Matrix<double> SquareMu::set_ab(unsigned int const& ref3) const {
	Matrix<double> tmp(2,2);
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
	return tmp;
}

unsigned int SquareMu::unit_cell_index(Vector<double> const& x) const {
	Vector<double> match(2,0);
	if(ref_(3)==2){ 
		if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 0; }
		match(0) = 0.2;
		match(1) = 0.4;
		if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 1; }
		match(0) = 0.4;
		match(1) = 0.8;
		if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 2; }
		match(0) = 0.6;
		match(1) = 0.2;
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
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return 5;
}
/*}*/

/*{method needed for checking*/
void SquareMu::lattice(){
	compute_H(0);

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

	double t;
	double mu;
	unsigned int s0;
	unsigned int s1;
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		xy0 = x_[s0];

		s1 = obs_[0](i,1);
		xy1 = x_[s1];

		mu = H_(s0,s0);
		if(std::abs(mu)>1e-4){
			if(mu<0){ color = "green"; }
			else    { color = "cyan"; }
			ps.circle(xy0,std::abs(mu),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
		}

		t = H_(s0,s1);
		if(obs_.size()>1){
			if(my::in_polygon(uc.row(),uc.ptr(),uc.ptr()+uc.row(),xy0(0),xy0(1))){ 
				t = obs_[1][obs_[0](i,2)].get_x(); 
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

			if(t>0){ color = "blue"; }
			else   { color = "red"; }

			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}

		if(obs_[0](i,3)){ ps.put(xy0(0)+0.1,(xy0(1)+xy1(1))/2.0,"\\tiny{"+std::string(1,my::int_to_alphabet(obs_[0](i,2),true))+"}"); }
		else            { ps.put((xy0(0)+xy1(0))/2.0,xy0(1)+0.1,"\\tiny{"+std::string(1,my::int_to_alphabet(obs_[0](i,2),true))+"}"); }

		if(i%2){ ps.put(xy0(0)+0.2,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}"); }
	}
	if(obs_.size()==4){
		double corr;
		double rescale(std::abs(0.25/obs_[3][1].get_x()));
		ps.cross(x_[0],0.25,"linecolor=black"); 
		ps.circle(x_[0],0.25,"linecolor=black"); 
		for(unsigned int i(1);i<n_;i++){
			corr = obs_[3][i].get_x();
			if(std::abs(corr)>1e-4){
				if(corr>0){ color = "blue"; }
				else      { color = "red"; }
				ps.circle(x_[i],sqrt(std::abs(corr*rescale)),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
			}
		}
	}
	ps.end(true,true,true);
}

void SquareMu::display_results(){
	lattice();

	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title("mu=" + my::tostring(mu_));
		std::string run_cmd("./mc -s:wf square-dimerizedbar");
		run_cmd += " -u:N " + my::tostring(N_);
		run_cmd += " -u:m " + my::tostring(m_);
		run_cmd += " -u:n " + my::tostring(n_);
		run_cmd += " -i:bc "+ my::tostring(bc_);
		run_cmd += " -d:mu "+ my::tostring(mu_);
		run_cmd += " -d -u:tmax 10";

		rst_file_->title(title,'-');
		rst_file_->change_text_onclick("run command",run_cmd);

		rst_file_->figure(dir_+filename_+".png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(dir_+filename_+".pdf")+RST::scale("200"));
	}
}

void SquareMu::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="square-mu";
	display_results();

	//plot_band_structure();
}
/*}*/
