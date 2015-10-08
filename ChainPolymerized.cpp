#include "ChainPolymerized.hpp"

ChainPolymerized::ChainPolymerized(System const& s, Vector<double> const& t):
	System(s),
	Chain<double>(set_spuc(t,N_/m_),"chain-polymerized"),
	t_(t)
{
	if(status_==2){
		init_fermionic();
		filename_ += "-t";
		for(unsigned int j(0);j<t_.size();j++){
			filename_ += ((t_(j)>0)?"+":"")+my::tostring(t_(j));
		}
		if(spuc_ != 1){
			system_info_.text("Trial wavefunction with different real hopping ");
			system_info_.text("terms for a "+RST::math("\\mathrm{SU}(N)")+" chain :"+RST::nl_);
			if(spuc_ != 4){
				std::string tmp("");
				for(unsigned int i(0);i<spuc_-1;i++){ tmp += "=\\bullet"; }
				system_info_.text(" "+RST::math("t_i : "+tmp+"-"));
			}
			else{system_info_.text(" "+RST::math("t_i : \\equiv\\bullet=\\bullet\\equiv\\bullet-"));}
			system_info_.nl();
		} else {
			system_info_.text("Trial wavefunction with uniform real hopping ");
			system_info_.text("terms for a "+RST::math("\\mathrm{SU}(N)")+" chain :");
		}
	}
}

/*{method needed for running*/
void ChainPolymerized::compute_H(){
	H_.set(n_,n_,0);
	Matrix<int> nb;
	for(unsigned int i(0);i<n_;i+=spuc_){
		for(unsigned int j(0);j<spuc_;j++){
			nb = get_neighbourg(i+j);
			H_(i+j,nb(0,0)) = nb(0,1)*t_(j);
		}
	}
	H_ += H_.transpose();
}

void ChainPolymerized::create(unsigned int const& which_observables){
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
	if(status_==2 && spuc_!=1){
		/*!Use the eigenvector (k1+k2)/sqrt(2) which correspond to the
		 * impulsion k1+k2=0.*/
		compute_H();
		diagonalize(false);
		double n1(0);
		double n2(0);
		std::complex<double> tmp1;
		std::complex<double> tmp2;
		unsigned int m(M_(0));
		for(unsigned int i(0);i<n_;i++){
			tmp1 = evec_(i,m) + evec_(i,m-1);//k=k1+k2=0
			tmp2 = evec_(i,m) - evec_(i,m-1);//k=k1-k2=2k1
			evec_(i,m-1)= tmp1;
			evec_(i,m)  = tmp2;
			n1 += my::norm_squared(tmp1);
			n2 += my::norm_squared(tmp2);
		}
		for(unsigned int i(0);i<n_;i++){
			evec_(i,m-1)/= sqrt(n1);
			evec_(i,m)  /= sqrt(n2);
		}
		for(unsigned int c(0);c<N_;c++){
			for(unsigned int i(0);i<n_;i++){
				EVec_[c](i,M_(c)-1) = real(evec_(i,M_(c)-1));
			}
		}
	}
}

void ChainPolymerized::save_param(IOFiles& w) const {
	std::string t_string("");
	for(unsigned int i(0);i<t_.size()-1;i++){
		t_string += my::tostring(t_(i))+",";
	}
	t_string += my::tostring(t_.back());
	w.write("t ("+t_string+")",t_);
}

unsigned int ChainPolymerized::set_spuc(Vector<double> const& t, unsigned int const& spuc){
	if(t.size() == spuc && !my::are_equal(t,Vector<double>(spuc,1.0))){ 
		return spuc; 
	} else { 
		std::cerr<<__PRETTY_FUNCTION__<<" : invalid t size : "<<t.size()<<std::endl;
		return 1; 
	}
}
/*}*/

/*{method needed for checking*/
void ChainPolymerized::check(){
	compute_H();
	plot_band_structure();
}
/*}*/

/*{method needed for analysing*/
std::string ChainPolymerized::extract_level_8(){
	rst_file_ = new RSTFile(info_+path_+dir_,filename_);
	std::string basename("../../../../../../../../"+analyse_+path_+dir_+filename_);
	std::string t_string("(");
	for(unsigned int i(0);i<t_.size()-1;i++){
		t_string += my::tostring(t_(i))+",";
	}
	t_string += my::tostring(t_.back())+")";
	std::string title("$N="+my::tostring(N_)+"$ $m="+my::tostring(m_)+"$ $n="+my::tostring(n_)+"$ bc="+my::tostring(bc_)+" $t_{ij}="+t_string+"$");

	/*!extract jdbin*/
	/*{*/
	IOFiles corr_file(analyse_+path_+dir_+filename_+"-corr.dat",true);
	IOFiles lr_corr_file(analyse_+path_+dir_+filename_+"-long-range-corr.dat",true);

	Vector<double> lr_corr_v(lr_corr_.size());
	Vector<double> poly_e(N_/m_,0);

	corr_file<<"%(2i+1)/2 corr(i,i+1) dx conv(0|1) #conv mean(0|1)"<<IOFiles::endl;
	lr_corr_file<<"%j corr(i,j) dx conv(0|1) #conv mean(0|1)"<<IOFiles::endl;

	(*data_write_)<<t_<<" "<<E_<<IOFiles::endl;
	for(unsigned int i(0);i<corr_.size();i++){
		corr_file<<i+0.5<<" "<<corr_[i]<<IOFiles::endl;
		poly_e(i%(N_/m_)) += corr_[i].get_x(); 
	}
	for(unsigned int i(0);i<lr_corr_.size();i++){
		lr_corr_file<<i<<" "<<lr_corr_[i]<<IOFiles::endl;
		lr_corr_v(i) = lr_corr_[i].get_x();
	}
	poly_e /= n_*m_/N_;
	poly_e.sort(std::less<double>());
	/*}*/
	/*!nearest neighbourg correlations*/
	/*{*/
	Gnuplot gp(analyse_+path_+dir_,filename_+"-corr");
	gp+="set key center";
	gp.label("x","site","offset 0,0.5");
	gp.label("y2","$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(i+1)>$");
	gp.title(title);
	gp+="plot '"+filename_+"-corr.dat' u 1:2:3 w errorbars lt 1 lc 7 notitle,\\";
	gp+="     "+my::tostring(poly_e(N_/m_-1)) + " w l lc 3 t 'd-merization="+my::tostring(poly_e(N_/m_-1)-poly_e(N_/m_-2))+"',\\";
	gp+="     "+my::tostring(poly_e(N_/m_-2)) + " w l lc 3 notitle";
	gp.save_file();
	//gp.create_image(true);
	rst_file_->figure(basename+"-corr.png","Correlation on links",RST::target(basename+"-corr.gp")+RST::width("1000"));
	/*}*/
	/*!long range correlations*/
	/*{*/
	unsigned int xi;
	unsigned int xf;
	Vector<double> exponents;
	bool fit(compute_critical_exponents(lr_corr_v,xi,xf,exponents));

	Gnuplot gplr(analyse_+path_+dir_,filename_+"-long-range-corr");
	gplr.range("x",N_/m_,n_-N_/m_);
	gplr.label("x","$\\|i-j\\|$","offset 0,0.5");
	gplr.label("y2","$<S_{\\alpha}^{\\alpha}(i)S_{\\alpha}^{\\alpha}(j)>-\\dfrac{m^2}{N}$","offset 1");
	gp.title(title);
	gplr+="set key center bottom";
	gplr+="set sample 1000";
	gplr+="m="+my::tostring(m_)+".0";
	gplr+="N="+my::tostring(N_)+".0";
	gplr+="n="+my::tostring(n_)+".0";
	gplr+="p0 = 1.0";
	gplr+="p1 = 2.0-2.0/N";
	gplr+="p2 = -1.0";
	gplr+="p3 = 2.0";
	gplr+="f(x) = p0*cos(2.0*pi*x*m/N)*(x**(-p1)+(n-x)**(-p1))+p2*(x**(-p3)+(n-x)**(-p3))";
	gplr+="set fit quiet";
	gplr+="fit [" + my::tostring(xi) + ":" + my::tostring(xf) + "] f(x) '"+filename_+"-long-range-corr.dat' u 1:2 noerrors via p0,p1,p2,p3"; 
	gplr+="plot '"+filename_+"-long-range-corr.dat' u 1:2:3 w errorbars lt 1 lc 7 notitle,\\";
	gplr+="     f(x) lc 7 " + std::string(fit?"lw 0.5":"dt 2") + " t sprintf('$\\eta=%f$, $\\mu=%f$',p1,p3)";
	gplr.save_file();
	//gplr.create_image(true);
	rst_file_->figure(basename+"-long-range-corr.png","Long range correlation",RST::target(basename+"-long-range-corr.gp")+RST::width("1000"));
	/*}*/
	/*!structure factor*/
	/*{*/
	unsigned int llr(lr_corr_.size());
	Vector<std::complex<double> > Ck(llr,0.0);
	std::complex<double> normalize(0.0);
	double dk(2.0*M_PI/llr);

	for(unsigned int k(0);k<llr;k++){
		for(unsigned int i(0);i<llr;i++){
			Ck(k) += std::polar(lr_corr_v(i),dk*k*i);
		}
		normalize += Ck(k); 
	}
	Ck /= dk*normalize;

	IOFiles data_sf(analyse_+path_+dir_+filename_+"-structure-factor.dat",true);
	for(unsigned int k(0);k<llr;k++){
		data_sf<<dk*k<<" "<<Ck(k).real()<<" "<<Ck(k).imag()<<IOFiles::endl;
	}

	Gnuplot gpsf(analyse_+path_+dir_,filename_+"-structure-factor");
	gp.title(title);
	gpsf+="set key bottom";
	gpsf.range("x","0","2*pi");
	switch(N_/m_){
		case 3: { gpsf+="set xtics ('0' 0,'$2\\pi/3$' 2.0*pi/3.0, '$4\\pi/3$' 4.0*pi/3.0,'$2\\pi$' 2.0*pi)"; } break;
		case 5: { gpsf+="set xtics ('0' 0,'$2\\pi/5$' 2.0*pi/5.0, '$4\\pi/5$' 4.0*pi/5.0, '$6\\pi/5$' 6.0*pi/5.0, '$8\\pi/5$' 8.0*pi/5.0, '$2\\pi$' 2.0*pi)"; } break;
		default:{ gpsf+="set xtics ('0' 0,'$\\pi/2$' pi/2.0,'$\\pi$' pi,'$3\\pi/2$' 3.0*pi/2.0,'$2\\pi$' 2.0*pi)"; } break;
	}
	gpsf.label("x","$k$","offset 0,0.5");
	gpsf.label("y2","$<S(k)>$");
	gpsf+="plot '"+filename_+"-structure-factor.dat' u 1:2 lt 1 lc 6 t 'real',\\";
	gpsf+="     '"+filename_+"-structure-factor.dat' u 1:3 lt 1 lc 7 t 'imag'";
	gpsf.save_file();
	//gpsf.create_image(true);
	rst_file_->figure(basename+"-structure-factor.png","Structure factor",RST::target(basename+"-structure-factor.gp")+RST::width("1000"));
	/*}*/
	/*!save*/
	/*{*/
	jd_write_->add_header()->title("System's parameters",'-');
	save_param(*jd_write_);
	save_input(*jd_write_);
	save_output(*jd_write_);
	jd_write_->write("polymerization strength",poly_e(N_/m_-1)-poly_e(N_/m_-2));
	jd_write_->write("critical exponents",exponents);
	/*}*/

	rst_file_->text(read_->get_header());
	rst_file_->save(false,true);
	delete rst_file_;
	rst_file_ = NULL;

	return t_string;
}

std::string ChainPolymerized::extract_level_7(){
	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp.title("$N="+my::tostring(N_)+"$ $m="+my::tostring(m_)+"$ $n="+my::tostring(n_)+"$ bc=$"+my::tostring(bc_)+"$");
	if(N_/m_!=4){
		gp.label("x","$t_"+my::tostring(N_/m_)+"$"," offset 0,1");
		gp.label("y2","$\\dfrac{E}{n}$","rotate by 0");
		gp.range("x","0.0","");
		gp+="f(x) = "+std::string(my::are_equal(t_(N_/m_-1),0)?"a+b*x**c":"a+b*(x-c)*(x-c)"); 
		gp+="a="+my::tostring(E_.get_x());
		gp+="b=1";
		gp+="c=1";
		gp+="set fit quiet";
		gp+="fit f(x) '"+filename_+".dat' u "+my::tostring(N_/m_)+":($" +my::tostring(N_/m_+5)+"==0?$"                              +my::tostring(N_/m_+1)+":1/0):"+my::tostring(N_/m_+2)+" zerror via a,b,c";
		gp+="plot '"+filename_+".dat' u "    +my::tostring(N_/m_)+":($" +my::tostring(N_/m_+5)+"==0?$"                              +my::tostring(N_/m_+1)+":1/0):"+my::tostring(N_/m_+2)+" lc 7 w e notitle";
		gp+="     f(x) lc 7 lw 0.5 "+std::string(my::are_equal(t_(N_/m_-1),0)?"notitle":"t sprintf('min %3.4f',c)");
	} else {
		gp.label("x","$t_2$");
		gp.label("y","$t_4$");
		gp.label("z","$\\dfrac{E}{n}$","rotate by 0");
		gp+="set xyplane 0";
		gp+="splot  '"+filename_+".dat' u 2:4:($9==0?$5:1/0) lc 7 notitle";
	}

	gp.save_file();
	gp.create_image(true,true);

	jd_write_->add_header()->title("System's parameters",'-');
	save_param(*jd_write_);
	save_input(*jd_write_);
	save_output(*jd_write_);
	jd_write_->write("polymerization strength",read_->read<double>());
	jd_write_->write("critical exponents",read_->read<Vector<double> >());

	return filename_;
}
/*}*/
