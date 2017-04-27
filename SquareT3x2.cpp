#include "SquareT3x2.hpp"

SquareT3x2::SquareT3x2(System const& s, Vector<double> const& t):
	System(s),
	Square<double>(set_ab(),6,"square-T3x2"),
	t_(t)
{
	if(t_.size()==12){
		if(status_==3){ init_lattice(); }
		if(status_==2){
			init_fermionic();

			system_info_.text("SquareT3x2 :");
			system_info_.item("Each color has the same Hamiltonian.");
			system_info_.item("6 sites in a 3x2 unit cell");

			filename_ += "-t";
			for(unsigned int i(0);i<t_.size();i++){
				filename_ += ((t_(i)>=0)?"+":"")+my::tostring(t_(i));
			}
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : t must contain 12 values (currently contains "<<t_.size()<<")"<<std::endl; }
}

/*{method needed for running*/
void SquareT3x2::compute_H(){
	H_.set(n_,n_,0);

	double t(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		switch(obs_[0](i,5)){
			case 0: { t = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)?1:0); }break;
			case 1: { t = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)?3:2); }break;
			case 2: { t = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)?5:4); }break;
			case 3: { t = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)?7:6); }break;
			case 4: { t = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)?9:8); }break;
			case 5: { t = (obs_[0](i,4)?bc_:1)*t_(obs_[0](i,3)?11:10); }break;

		}
		H_(obs_[0](i,0),obs_[0](i,1)) = t;
	}
	H_ += H_.transpose();
}

void SquareT3x2::create(){
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

void SquareT3x2::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("t=(");
		Vector<double> param(t_.size());

		for(unsigned int i(0);i<t_.size()-1;i++){
			param(i) = t_(i);
			s += my::tostring(t_(i))+",";
		}
		param(t_.size()-1) = t_.back();
		s += my::tostring(t_.back())+")";

		w.add_to_header()->title(s,'<');
		w<<param;
		w.add_to_header()->add(system_info_.get());
	} else { w<<t_<<" "; }
}

Matrix<double> SquareT3x2::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 3.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 0.0;
	tmp(1,1) = 2.0;
	return tmp;
}

unsigned int SquareT3x2::unit_cell_index(Vector<double> const& x) const {
	if(my::are_equal(x(1),0.0,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.0,    eq_prec_,eq_prec_)){ return 0; }
		if(my::are_equal(x(0),1.0/3.0,eq_prec_,eq_prec_)){ return 1; }
		if(my::are_equal(x(0),2.0/3.0,eq_prec_,eq_prec_)){ return 2; }
	}
	if(my::are_equal(x(1),0.5,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.0,    eq_prec_,eq_prec_)){ return 3; }
		if(my::are_equal(x(0),1.0/3.0,eq_prec_,eq_prec_)){ return 4; }
		if(my::are_equal(x(0),2.0/3.0,eq_prec_,eq_prec_)){ return 5; }
	}
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void SquareT3x2::display_results(){
	compute_H();
	std::string t(my::tostring(t_));
	draw_lattice(true,true,false,(dir_nn_[2]+dir_nn_[3])*0.5,"-d:t "+t,"t=("+t+")");
}

void SquareT3x2::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="square-T3x2";
	display_results();

	//compute_H();
	//plot_band_structure();
}

void SquareT3x2::param_fit_therm_limit(std::string& f, std::string& param, std::string& range){
	f="f(x) = a+b*x*x";
	param = "a,b";
	range = "[0:0.015]";
}

std::string SquareT3x2::extract_level_2(){
	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp+="f(x) = a+b*x*x";
	gp+="g(x,a,b) = a+b*x*x";
	gp+="set fit quiet";
	gp+="fit f(x) '"+filename_+".dat' u (1.0/$3):($5/($1*$1)):($6/($1*$1)) every 2::1 yerror via a,b";
	gp+="a0=a";
	gp+="b0=b";
	gp+="fit f(x) '"+filename_+".dat' u (1.0/$3):($5/($1*$1)):($6/($1*$1)) every 2::2 yerror via a,b";
	gp+="a1=a";
	gp+="b1=b";
	gp+="set print \'../"+sim_.substr(0,sim_.size()-1)+".dat\' append";
	gp+="print " + my::tostring(N_) + "," + my::tostring(m_) + ",a0,a1,(a0+a1)/2.0,mymin(a0,a1)";
	gp.range("x","0","");
	gp.label("x","$\\frac{ 1}{n}$");
	gp.label("y2","$\\frac{E}{nN^2}$","rotate by 0");
	gp+="plot '"+filename_+".dat' u (1.0/$3):($5/($1*$1)):($6/($1*$1)) w e lc 1 notitle,\\";
	gp+="     '"+filename_+".dat' u (1.0/$3):($5/($1*$1)):($6/($1*$1)) every 2::1 w e lc 2 notitle,\\";
	gp+="     '"+filename_+".dat' u (1.0/$3):($5/($1*$1)):($6/($1*$1)) every 2::2 w e lc 3 notitle,\\";
	gp+="     [0:0.015] g(x,a0,b0) lc 2 t sprintf('$a_0=%f$',a0),\\";
	gp+="     [0:0.015] g(x,a1,b1) lc 3 t sprintf('$a_1=%f$',a1)";
	gp.save_file();
	gp.create_image(true,"png");

	return filename_;
}
/*}*/
