#include"LadderFreeFlux.hpp"

LadderFreeFlux::LadderFreeFlux(System const& s, Vector<double> const& t, Vector<double> const& flux):
	System(s),
	Ladder<std::complex<double> >(set_spuc(t,flux),"ladder-freeflux"),
	t_(t),
	flux_(flux)
{
	if(status_==2){
		init_fermionic();

		system_info_.text("LadderFreeFlux :");
		system_info_.item("Each color has the same Hamiltonian.");
		system_info_.item("Free hopping term");
		system_info_.item("Free flux per plaquette ");
		filename_ += "-t";
		for(unsigned int i(0);i<t_.size();i++){
			filename_ += ((t_(i)>0)?"+":"")+my::tostring(t_(i));
		}
		filename_ += "-phi";
		for(unsigned int i(0);i<flux_.size()-1;i++){
			filename_ += ((flux_(i)>0)?"+":"")+my::tostring(flux_(i)); 
		}
	}
}

/*{method needed for running*/
void LadderFreeFlux::compute_H(){
	H_.set(n_,n_,0);
	Matrix<int> nb;
	std::cerr<<__PRETTY_FUNCTION__<<" : undefined"<<std::endl;
	//for(unsigned int i(0);i<n_;i++){
		//nb = get_neighbourg(i);
		//if(i%spuc_){
			//if(i%2){
				//H_(i,nb(0,0)) = my::chop(std::polar(nb(0,1)*t_(k++),flux_(f++)*M_PI));
			//} else {
				//H_(i,nb(0,0)) = nb(0,1)*t_(k++);
				//H_(i,nb(1,0)) = nb(1,1)*t_(k++);
			//}
		//} else {
			///*!decoupled chains in the limit J⊥(0)=J_(1)->0, therefore, in
			// * that case the variational parameter is t⊥ (t‖ is set to 1),
			// * otherwise the inverse is done*/
			//if(J_(0)>J_(1)){
				//H_(i,nb(0,0)) = nb(0,1);
				//H_(i,nb(1,0)) = nb(1,1)*t_(k++);
			//} else {
				//H_(i,nb(0,0)) = nb(0,1)*t_(k++);
				//H_(i,nb(1,0)) = nb(1,1);
			//}
		//}
		//k = k%t_.size();
		//f = f%flux_.size();
	//}
	H_ += H_.conjugate_transpose();
}

void LadderFreeFlux::create(){
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

void LadderFreeFlux::save_param(IOFiles& w) const {
	if(w.is_binary()){
		Vector<double> param(t_.size()+flux_.size());
		std::string s("t=(");
		for(unsigned int i(0);i<t_.size()-1;i++){ 
			s += my::tostring(t_(i))+","; 
			param(i) = t_(i);
		}
		param(t_.size()-1) = t_.back();
		s += my::tostring(t_.back())+") : " + RST::math("\\phi/\\pi")+"=(";
		for(unsigned int i(0);i<flux_.size()-1;i++){
			s += my::tostring(flux_(i))+","; 
			param(i+t_.size()) = flux_(i);
		}
		param.back() = flux_.back();
		s += my::tostring(flux_.back())+")";
		w.add_header()->title(s,'<');
		w<<param;
		GenericSystem<std::complex<double> >::save_param(w);
	} else { w<<t_<<" "<<flux_<<" "; }
}

unsigned int LadderFreeFlux::set_spuc(Vector<double> const& t, Vector<double> const& flux){
	if((t.size()+1)%flux.size()){
		std::cerr<<__PRETTY_FUNCTION__<<" : incoherent t and flux sizes : "<<t.size()<<" "<<flux.size()<<std::endl;
		return 1;
	} else {
		switch(t.size()){
			case 2: { return 2; } break;
			case 5: { return 4; } break;
			case 8: { return 6; } break;
			case 11:{ return 8; } break;
			default:{
						std::cerr<<__PRETTY_FUNCTION__<<" : invalid t size : "<<t.size()<<std::endl;
						return 1;
					}
		}
	}
}

void LadderFreeFlux::get_wf_symmetries(std::vector<Matrix<int> >& sym) const {
	int A(J_(0)>J_(1)?0:-1); //link between sites 0-1 (rung)
	int B(J_(0)>J_(1)?-1:0); //link between sites 0-2 (leg)
	switch(spuc_){
		case 8:
			{
				Matrix<int> tmp;
				/*{no symmetry breaking*/
				sym.push_back(tmp);
				/*}*/

				/*{facing tetramerization*/
				tmp.set(7,3);
				tmp(0,0) = 1;
				tmp(0,1) = B;
				tmp(0,2) = 1;

				tmp(1,0) = 4;
				tmp(1,1) = 2;
				tmp(1,2) = 1;

				tmp(2,0) = 5;
				tmp(2,1) = B;
				tmp(2,2) = 1;

				tmp(3,0) = 6;
				tmp(3,1) = 3;
				tmp(3,2) = 1;

				tmp(4,0) = 7;
				tmp(4,1) = B;
				tmp(4,2) = 1;

				tmp(5,0) = 9;
				tmp(5,1) = A;
				tmp(5,2) = 1;

				tmp(6,0) = 10;
				tmp(6,1) = 8;
				tmp(6,2) = 1;
				sym.push_back(tmp);
				/*}*/

				/*{shifted by one tetramerization*/
				tmp.set(6,3);
				tmp(0,0) = 1;
				tmp(0,1) = 2;
				tmp(0,2) = 1;

				tmp(1,0) = 4;
				tmp(1,1) = 5;
				tmp(1,2) = 1;

				tmp(2,0) = 7;
				tmp(2,1) = 2;
				tmp(2,2) = 1;

				tmp(3,0) = 8;
				tmp(3,1) = 2;
				tmp(3,2) = 1;

				tmp(4,0) = 9;
				tmp(4,1) = 3;
				tmp(4,2) = 1;

				tmp(5,0) = 10;
				tmp(5,1) = B;
				tmp(5,2) = 1;
				sym.push_back(tmp);
				/*}*/

				/*{shifted by two tetramerization*/
				tmp.set(8,3);
				tmp(0,0) = 1;
				tmp(0,1) = B;
				tmp(0,2) = 1;

				tmp(1,0) = 3;
				tmp(1,1) = A;
				tmp(1,2) = 1;

				tmp(2,0) = 4;
				tmp(2,1) = 8;
				tmp(2,2) = 1;

				tmp(3,0) = 5;
				tmp(3,1) = B;
				tmp(3,2) = 1;

				tmp(4,0) = 6;
				tmp(4,1) = A;
				tmp(4,2) = 1;

				tmp(5,0) = 7;
				tmp(5,1) = B;
				tmp(5,2) = 1;

				tmp(6,0) = 9;
				tmp(6,1) = A;
				tmp(6,2) = 1;

				tmp(7,0) = 10;
				tmp(7,1) = 2;
				tmp(7,2) = 1;
				sym.push_back(tmp);
				/*}*/

				/*{double dimerization*/
				tmp.set(7,3);
				tmp(0,0) = 1;
				tmp(0,1) = B;
				tmp(0,2) = 1;

				tmp(1,0) = 2;
				tmp(1,1) = B;
				tmp(1,2) = 1;

				tmp(2,0) = 4;
				tmp(2,1) = B;
				tmp(2,2) = 1;

				tmp(3,0) = 6;
				tmp(3,1) = A;
				tmp(3,2) = 1;

				tmp(4,0) = 7;
				tmp(4,1) = 5;
				tmp(4,2) = 1;

				tmp(5,0) = 8;
				tmp(5,1) = 5;
				tmp(5,2) = 1;

				tmp(6,0) = 10;
				tmp(6,1) = 5;
				tmp(6,2) = 1;
				sym.push_back(tmp);
				/*}*/
			}break;
		default:{ std::cerr<<__PRETTY_FUNCTION__<<"unknown spuc_"<<std::endl; }
	}
}
/*}*/

/*{method needed for checking*/
void LadderFreeFlux::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="ladder-freeflux";
	display_results();

	//plot_band_structure();
}

void LadderFreeFlux::plot(bool const& create_image){
	if(obs_.size()>3){
		/*!long range correlations*/
		/*{*/
		unsigned int llr(obs_[1].nval());
		Vector<std::complex<double> > Ck_intra(llr,0.0);
		Vector<std::complex<double> > Ck_inter(llr,0.0);
		std::complex<double> normalize_intra(0.0);
		std::complex<double> normalize_inter(0.0);
		double dk(2.0*M_PI/llr);

		for(unsigned int k(0);k<llr;k++){
			for(unsigned int i(0);i<llr;i++){
				Ck_intra(k) += std::polar(obs_[1][i].get_x(),dk*k*i);
				Ck_inter(k) += std::polar(obs_[2][i].get_x(),dk*k*i);
			}
			normalize_intra += Ck_intra(k);
			normalize_inter += Ck_inter(k);
		}
		Ck_intra /= dk*normalize_intra;
		Ck_inter /= dk*normalize_inter;

		IOFiles file_c(analyse_+path_+dir_+filename_+"-lr-c.dat",true);
		IOFiles file_sf(analyse_+path_+dir_+filename_+"-lr-sf.dat",true);
		for(unsigned int l(0);l<llr;l++){
			file_c<<l<<" "<<obs_[1][l]<<" "<<obs_[2][l]<<IOFiles::endl;
			file_sf<<dk*l<<" "<<Ck_intra(l).real()<<" "<<Ck_intra(l).imag()<<" "<<Ck_inter(l).real()<<" "<<Ck_inter(l).imag()<<IOFiles::endl;
		}

		Gnuplot gp(analyse_+path_+dir_,filename_+"-lr");
		gp.multiplot();
		/*{correlations*/
		gp.range("x","0",llr/2);

		gp.tics("x");
		gp.margin("0.1","0.5","0.9","0.5");
		gp+="plot '"+filename_+"-lr-c.dat' u 1:2:3 w errorbars lt 1 lc 6 t 'intra'";

		gp.margin("0.1","0.5","0.5","0.1");
		gp.tics("x","");
		gp+="plot '"+filename_+"-lr-c.dat' u 1:6:7  w errorbars lt 1 lc 7 t 'inter'";
		/*}*/
		/*{structure factor*/
		gp.range("x","0","pi");
		//gp.range("y2","0","");

		gp.key("left");
		gp.tics("x");
		gp.tics("y");
		gp.tics("x2","('' pi/3, '' pi/2, '' 2*pi/3, '' pi) mirror");
		gp.tics("y2","mirror");
		gp.margin("0.5","0.9","0.9","0.5");
		gp+="plot '"+filename_+"-lr-sf.dat' u 1:2 axes x1y2 lt 1 lc 6 notitle,\\";
		gp+="     '"+filename_+"-lr-sf.dat' u 1:3 axes x1y2 lt 2 lc 6 notitle";

		gp.margin("0.5","0.9","0.5","0.1");
		gp.tics("x","('$\\pi/3$' pi/3, '$\\pi/2$' pi/2, '$2\\pi/3$' 2*pi/3, '$\\pi$' pi)");
		gp+="plot '"+filename_+"-lr-sf.dat' u 1:4 axes x1y2 lt 1 lc 7 notitle,\\";
		gp+="     '"+filename_+"-lr-sf.dat' u 1:5 axes x1y2 lt 2 lc 7 notitle";
		/*}*/
		gp.save_file();
		gp.create_image(create_image,true);
		/*}*/
	}
	if(obs_.size()==5){
		/*!(anti)symmetric correlations and structure factors*/
		/*{*/
		unsigned int llr(obs_[2].nval());
		Vector<std::complex<double> > Ckm(llr,0.0);
		Vector<std::complex<double> > Ckp(llr,0.0);
		std::complex<double> normalize_m(0.0);
		std::complex<double> normalize_p(0.0);
		double dk(2.0*M_PI/llr);

		for(unsigned int k(0);k<llr;k++){
			for(unsigned int i(0);i<llr;i++){
				Ckm(k) += std::polar(obs_[1][i].get_x()-obs_[2][i].get_x()-obs_[3][i].get_x()+obs_[4][i].get_x(),dk*k*i);
				Ckp(k) += std::polar(obs_[1][i].get_x()+obs_[2][i].get_x()+obs_[3][i].get_x()+obs_[4][i].get_x(),dk*k*i);
			}
			normalize_m += Ckm(k);
			normalize_p += Ckp(k);
		}
		Ckm /= dk*normalize_m;
		Ckp /= dk*normalize_p;

		IOFiles file_c(analyse_+path_+dir_+filename_+"-as-c.dat",true);
		IOFiles file_sf(analyse_+path_+dir_+filename_+"-as-sf.dat",true);
		for(unsigned int l(0);l<llr;l++){
			file_c<<l<<" "<<obs_[1][l]<<" "<<obs_[2][l]<<" "<<obs_[3][l]<<" "<<obs_[4][l]<<IOFiles::endl;
			file_sf<<dk*l<<" "<<Ckm(l).real()<<" "<<Ckm(l).imag()<<" "<<Ckp(l).real()<<" "<<Ckp(l).imag()<<IOFiles::endl;
		}

		Gnuplot gp(analyse_+path_+dir_,filename_+"-as");
		gp.multiplot();
		/*{correlations*/
		gp.range("x","0",llr/2);

		gp.tics("x");
		gp.margin("0.1","0.5","0.9","0.5");
		gp+="plot '"+filename_+"-as-c.dat' u 1:($2-$6-$10+$14):($3+$7+$18+$18) w errorbars lt 1 lc 6 t '$-$'";

		gp.margin("0.1","0.5","0.5","0.1");
		gp.tics("x","");
		gp+="plot '"+filename_+"-as-c.dat' u 1:($2+$6+$10+$14):($3+$7+$18+$18)  w errorbars lt 1 lc 7 t '$+$'";
		/*}*/
		/*{structure factor*/
		gp.range("x","0","pi");
		//gp.range("y2","0","");

		gp.key("left");
		gp.tics("x");
		gp.tics("y");
		gp.tics("x2","('' pi/3, '' pi/2, '' 2*pi/3, '' pi)");
		gp.tics("y2","");
		gp.margin("0.5","0.9","0.9","0.5");
		gp+="plot '"+filename_+"-as-sf.dat' u 1:2 axes x1y2 lt 1 lc 6 notitle,\\";
		gp+="     '"+filename_+"-as-sf.dat' u 1:3 axes x1y2 lt 2 lc 6 notitle";

		gp.margin("0.5","0.9","0.5","0.1");
		gp.tics("x","('$\\pi/3$' pi/3, '$\\pi/2$' pi/2, '$2\\pi/3$' 2*pi/3, '$\\pi$' pi)");
		gp+="plot '"+filename_+"-as-sf.dat' u 1:4 axes x1y2 lt 1 lc 7 notitle,\\";
		gp+="     '"+filename_+"-as-sf.dat' u 1:5 axes x1y2 lt 2 lc 7 notitle";
		/*}*/
		gp.save_file();
		gp.create_image(create_image,true);
		/*}*/
	}
}

void LadderFreeFlux::lattice(){
	compute_H();
	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	double x_shift(spuc_/2+2);

	PSTricks ps(info_+path_+dir_,filename_+"-pstricks");
	ps.begin(-1,-5,n_/1.5,2,filename_+"-pstricks");
	std::complex<double> t;
	double corr;
	unsigned int s0;
	unsigned int s1;
	for(unsigned int i(0);i<3*spuc_/2;i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		xy0(0) = s0/2;
		xy0(1) = s0%2;
		xy1(0) = s1/2;
		xy1(1) = s1%2;

		t=H_(s0,s1);
		if(std::abs(t)>1e-4){
			if(xy1(0)<xy0(0)){
				xy1(0) = xy0(0)+1;
				linestyle="dashed";
			} else { linestyle="solid"; }

			if(t.real()>0){ color = "blue"; }
			else          { color = "red"; }
			linewidth = my::tostring(std::abs(t))+"mm";

			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1),"linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
			switch(i%3){
				case 0: { ps.put((xy0(0)+xy1(0))/2.0,xy0(1)+0.2,"\\tiny{"+my::tostring(t)+"}"); }break;
				case 1: { ps.put(xy0(0)+0.2,(xy0(1)+xy1(1))/2.0,"\\tiny{"+my::tostring(t)+"}"); }break;
				case 2: { ps.put((xy0(0)+xy1(0))/2.0,xy0(1)-0.2,"\\tiny{"+my::tostring(t)+"}"); }break;
			}
		}

		if(obs_[0].nval()){/*bound energy*/
			corr = obs_[0][i].get_x();
			if(std::abs(corr)>1e-4){
				if(corr<0){ color = "red"; }
				else { color = "blue"; }
				linewidth = my::tostring(std::abs(corr))+"mm";

				ps.line("-",xy0(0)+x_shift,xy0(1),xy1(0)+x_shift,xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
			}
			if(i%3==0){ 
				ps.put(xy0(0)+x_shift,xy0(1)-0.2,"\\tiny{"+my::tostring(s0)+"}"); 
				ps.put((xy0(0)+xy1(0))/2.0+2*x_shift,xy0(1),"\\tiny{"+my::tostring(corr).substr(0,5)+"}");
			}
			if(i%3==1){ 
				ps.put(xy1(0)+x_shift,xy1(1)+0.2,"\\tiny{"+my::tostring(s1)+"}"); 
				ps.put((xy0(0)+xy1(0))/2.0+2*x_shift,(xy0(1)+xy1(1))/2.0,"\\tiny{"+my::tostring(corr).substr(0,5)+"}");
			}
			if(i%3==2){
				ps.put((xy0(0)+xy1(0))/2.0+2*x_shift,xy1(1),"\\tiny{"+my::tostring(corr).substr(0,5)+"}");
			}
		}

		if(i%3==0){ ps.put(xy0(0),xy0(1)-0.2,"\\tiny{"+my::tostring(s0)+"}"); }
		if(i%3==1){ ps.put(xy1(0),xy1(1)+0.2,"\\tiny{"+my::tostring(s1)+"}"); }
	}
	if(obs_.size()==5){/*long range correlations*/
		double rescale(0.75/obs_[1][0].get_x());
		unsigned int n(std::min(4*N_/m_+1,obs_[1].nval()));
		unsigned int idx;
		for(unsigned int i(0);i<n;i++){
			idx = (obs_[1].nval()-n/2+i)%obs_[1].nval();

			corr = obs_[1][idx].get_x()*rescale;
			xy0(0) = i;
			xy0(1) = -2;
			if(std::abs(corr)>1e-4){
				if(i!=idx){
					if(corr<0){ color = "red"; }
					else { color = "blue"; }
				} else { color = "black"; }
				ps.circle(xy0,std::abs(corr),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
			}

			corr = obs_[2][idx].get_x()*rescale;
			xy0(1) = -1;
			if(std::abs(corr)>1e-4){
				if(corr<0){ color = "red"; }
				else { color = "blue"; }
				ps.circle(xy0,std::abs(corr),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
			}
		}
	}
	ps.end(true,true,true);
}

void LadderFreeFlux::display_results(){
	lattice();
	plot(true);
	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }
		std::string title(RST::math("\\theta=")+my::tostring(acos(this->J_(0))) + " : t=(");
		for(unsigned int i(0);i<t_.size()-1;i++){ title += my::tostring(t_(i)) + ","; }
		title += my::tostring(t_.back()) + ") "+RST::math("\\phi/\\pi=")+"(";
		for(unsigned int i(0);i<flux_.size()-1;i++){ title += my::tostring(flux_(i)) + ","; }
		title += my::tostring(flux_.back()) + ")";
		if(dir_ == "P/" || dir_ == "O/" || dir_ == "A/"){
			rst_file_->title("|theta"+my::tostring(acos(this->J_(0)))+"|_",'-');
			rst_file_->replace("theta"+my::tostring(acos(this->J_(0))),title);
		} else { rst_file_->title(title,'-'); }
		rst_file_->figure(dir_+filename_+"-pstricks.png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(dir_+filename_+"-pstricks.pdf")+RST::scale("200"));
		rst_file_->figure(relative_path+filename_+"-lr.png","long range correlations",RST::target(relative_path+filename_+"-lr.gp")+RST::scale("200"));
		rst_file_->figure(relative_path+filename_+"-as.png","(anti)symmetric correlations",RST::target(relative_path+filename_+"-as.gp")+RST::scale("200"));
	}
}
/*}*/

/*{method needed for analysing*/
std::string LadderFreeFlux::extract_level_6(){
	(*data_write_)<<N_<<" "<<m_<<" "<<n_<<" "<<bc_<<" "<<asin(J_(1))<<" "<<obs_[0][0]<<IOFiles::endl;

	display_results();

	save_param(*jd_write_);
	save(*jd_write_);

	return filename_;
}
/*}*/
