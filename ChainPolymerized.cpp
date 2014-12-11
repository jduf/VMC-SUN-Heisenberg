#include "ChainPolymerized.hpp"

ChainPolymerized::ChainPolymerized(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc, double delta):
	System(ref,N,m,n,M,bc),
	Chain<double>(are_equal(delta,0)?1:N_/m_,"chain-polymerized"),
	delta_(delta)
{
	if(status_==1){
		init_fermionic();
		filename_ += "-delta" + tostring(delta_);
		if(!are_equal(delta,0)){
			system_info_.text("+ Spin chain, with different real hopping term.");
			system_info_.text("  For N colors and m particules per sites, every");
			system_info_.text("  N/m, there is a weaker bound, namely t-delta");
			system_info_.text("  instead of t+delta. (t=1,delta>0) :");
		} else {
			system_info_.text("+ Spin chain with real and identical hopping");
			system_info_.text("  term between each sites :");
		}
	}
}

/*{method needed for running*/
void ChainPolymerized::compute_H(){
	/*!If t<0, delta<0 otherwise no polymerization occurs
	 * If t>0, delta>0 otherwise no polymerization occurs */
	double t(1.0);
	H_.set(n_,n_,0);
	Matrix<int> nb;
	unsigned int a(n_/L_);
	for(unsigned int i(0); i<n_; i+=a){
		for(unsigned int j(0); j<a-1; j++){
			nb = get_neighbourg(i+j);
			H_(i+j,nb(0,0)) = t+delta_;
		}
		nb = get_neighbourg(i+a-1);
		H_(i+a-1,nb(0,0)) = nb(0,1)*(t-delta_);
	}
	H_ += H_.transpose();
}

void ChainPolymerized::create(){
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);
	lr_corr_.set(links_.row(),50,5,false);

	compute_H();
	diagonalize_H(H_);
	for(unsigned int c(0);c<N_;c++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_(c);j++){
				EVec_[c](i,j) = H_(i,j);
			}
		}
	}
	if(degenerate_){
		degenerate_ = false;
		compute_H();
		select_eigenvectors(M_(0));
		for(unsigned int c(0);c<N_;c++){
			for(unsigned int i(0);i<n_;i++){
				EVec_[c](i,M_(c)-1) = real(evec_(i,M_(c)-1));
			}
		}
	}
}

void ChainPolymerized::save() const {
	GenericSystem<double>::save();
	jd_write_->write("delta (t+-delta)",delta_);
}
/*}*/

/*{method needed for checking*/
void ChainPolymerized::check(){
	compute_H();
	plot_band_structure();
}
/*}*/

/*{method needed for analysing*/
std::string ChainPolymerized::extract_level_7(){
	rst_file_ = new RSTFile(info_+path_+dir_,filename_);

	/*!extract jdbin*/
	/*{*/
	IOFiles corr_file(analyse_+path_+dir_+filename_+"-corr.dat",true);
	IOFiles lr_corr_file(analyse_+path_+dir_+filename_+"-long-range-corr.dat",true);

	Vector<double> lrc_mean(links_.row(),0);
	Vector<double> poly_e(N_/m_,0);
	unsigned int nruns;
	unsigned int tmax;

	(*read_)>>nruns>>tmax;
	(*data_write_)<<"%delta E dx conv(0|1) #conv mean(0|1)"<<IOFiles::endl;
	corr_file<<"%(2i+1)/2 corr(i,i+1) dx conv(0|1) #conv mean(0|1)"<<IOFiles::endl;
	lr_corr_file<<"%j corr(i,j) dx conv(0|1) #conv mean(0|1)"<<IOFiles::endl;
	/*!the +1 is the average over all runs */
	for(unsigned int i(0);i<nruns+1;i++){ 
		(*read_)>>E_>>corr_>>lr_corr_;
		(*data_write_)<<delta_<<" "<<E_<<" "<<(i<nruns)<<IOFiles::endl;
		for(unsigned int j(0);j<corr_.size();j++){
			corr_file<<j+0.5<<" "<<corr_[j]<<" "<<(i<nruns)<<IOFiles::endl;
		}
		for(unsigned int j(0);j<lr_corr_.size();j++){
			lr_corr_file<<j<<" "<<lr_corr_[j]<<" "<<(i<nruns)<<IOFiles::endl;
		}
		if(i<nruns){
			for(unsigned int j(0);j<lr_corr_.size();j++){
				lrc_mean(j) += lr_corr_[j].get_x()/nruns;
			}
			unsigned int k(0);
			do{ poly_e(k%(N_/m_)) += corr_[k].get_x(); }
			while(++k<corr_.size());
		} else {
			for(unsigned int j(0);j<lr_corr_.size();j++){
				if(lr_corr_[j].get_conv()){ lrc_mean(j) = lr_corr_[j].get_x(); } 
			}
		}
	}
	poly_e /= nruns*n_*m_/N_;
	poly_e.sort(std::less<double>());
	/*}*/
	/*!nearest neighbourg correlations*/
	/*{*/
	Gnuplot gp(analyse_+path_+dir_,filename_+"-corr");
	gp+="set xlabel 'site' offset 0,0.5";
	gp+="set y2label '$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(i+1)>$'";
	gp+="set key center";
	gp+="set title '$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$ bc="+tostring(bc_)+" $\\delta="+tostring(delta_)+"$'";
	gp+="plot '"+filename_+"-corr.dat' u 1:(($6==1 && $5==0)?$2:1/0):3 w errorbars lt 1 lc 5 t 'Not converged',\\";
	gp+="     '"+filename_+"-corr.dat' u 1:(($6==1 && $5==1)?$2:1/0):3 w errorbars lt 1 lc 6 t 'Converged',\\";
	gp+="     '"+filename_+"-corr.dat' u 1:($6==0?$2:1/0):3 w errorbars lt 1 lc 7 t 'Mean',\\";
	gp+="     "+tostring(poly_e(N_/m_-1)) + " w l lc 3 t 'd-merization="+tostring(poly_e(N_/m_-1)-poly_e(N_/m_-2))+"',\\";
	gp+="     "+tostring(poly_e(N_/m_-2)) + " w l lc 3 notitle";
	gp.save_file();
	//gp.create_image(true);
	rst_file_->link_figure(analyse_+path_+dir_+filename_+"-corr.png","Correlation on links",analyse_+path_+dir_+filename_+"-corr.gp",1000);
	/*}*/
	/*!long range correlations*/
	/*{*/
	unsigned int xi(1);
	unsigned int xf(n_);
	Vector<double> exponents(4);
	compute_critical_exponents(xi,xf,exponents,lrc_mean);

	Gnuplot gplr(analyse_+path_+dir_,filename_+"-long-range-corr");
	gplr.xrange(N_/m_,n_-N_/m_);
	gplr+="set xlabel '$\\|i-j\\|$' offset 0,0.5";
	gplr+="set y2label '$<S_{\\alpha}^{\\alpha}(i)S_{\\alpha}^{\\alpha}(j)>-\\dfrac{m^2}{N}$' offset 1";
	gplr+="set title '$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$ bc="+tostring(bc_)+" $\\delta="+tostring(delta_)+"$'";
	gplr+="set key center bottom";
	gplr+="set sample 1000";
	gplr+="m="+tostring(m_)+".0";
	gplr+="N="+tostring(N_)+".0";
	gplr+="n="+tostring(n_)+".0";
	gplr+="p0 = 1.0";
	gplr+="p1 = 2.0-2.0/N";
	gplr+="p2 = -1.0";
	gplr+="p3 = 2.0";
	gplr+="f(x) = p0*cos(2.0*pi*x*m/N)*(x**(-p1)+(n-x)**(-p1))+p2*(x**(-p3)+(n-x)**(-p3))";
	gplr+="set fit quiet";
	gplr+="fit [" + tostring(xi) + ":" + tostring(xf) + "] f(x) '"+filename_+"-long-range-corr.dat' u 1:($6==0?$2:1/0) noerrors via p0,p1,p2,p3"; 
	gplr+="plot '"+filename_+"-long-range-corr.dat' u 1:(($6==1 && $5==0)?$2:1/0):3 w errorbars lt 1 lc 5 t 'Not converged',\\";
	gplr+="     '"+filename_+"-long-range-corr.dat' u 1:(($6==1 && $5==1)?$2:1/0):3 w errorbars lt 1 lc 6 t 'Converged',\\";
	gplr+="     '"+filename_+"-long-range-corr.dat' u 1:($6==0?$2:1/0):3 w errorbars lt 1 lc 7 t 'Mean',\\";
	gplr+="     f(x) lc 7 lw 0.5 t sprintf('$\\eta=%f$, $\\mu=%f$',p1,p3)";
	gplr.save_file();
	//gplr.create_image(true);
	rst_file_->link_figure(analyse_+path_+dir_+filename_+"-long-range-corr.png","Long range correlation",analyse_+path_+dir_+filename_+"-long-range-corr.gp",1000);
	/*}*/
	/*!structure factor*/
	/*{*/
	unsigned int llr(lr_corr_.size());
	Vector<std::complex<double> > Ck(llr,0.0);
	std::complex<double> normalize(0.0);
	double dk(2.0*M_PI/llr);

	for(unsigned int k(0);k<llr;k++){
		for(unsigned int i(0);i<llr;i++){
			Ck(k) += std::polar(lrc_mean(i),dk*k*i);
		}
		normalize += Ck(k); 
	}
	Ck /= dk*normalize;

	IOFiles data_sf(analyse_+path_+dir_+filename_+"-structure-factor.dat",true);
	for(unsigned int k(0);k<llr;k++){
		data_sf<<dk*k<<" "<<Ck(k).real()<<" "<<Ck(k).imag()<<IOFiles::endl;
	}

	Gnuplot gpsf(analyse_+path_+dir_,filename_+"-structure-factor");
	gpsf+="set title '$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$ bc="+tostring(bc_)+" $\\delta="+tostring(delta_)+"$'";
	gpsf+="set key bottom";
	gpsf.xrange("0","2*pi");
	switch(N_/m_){
		case 3: { gpsf+="set xtics ('0' 0,'$2\\pi/3$' 2.0*pi/3.0, '$4\\pi/3$' 4.0*pi/3.0,'$2\\pi$' 2.0*pi)"; } break;
		case 5: { gpsf+="set xtics ('0' 0,'$2\\pi/5$' 2.0*pi/5.0, '$4\\pi/5$' 4.0*pi/5.0, '$6\\pi/5$' 6.0*pi/5.0, '$8\\pi/5$' 8.0*pi/5.0, '$2\\pi$' 2.0*pi)"; } break;
		default:{ gpsf+="set xtics ('0' 0,'$\\pi/2$' pi/2.0,'$\\pi$' pi,'$3\\pi/2$' 3.0*pi/2.0,'$2\\pi$' 2.0*pi)"; } break;
	}
	gpsf+="set xlabel '$k$' offset 0,0.5";
	gpsf+="set y2label '$<S(k)>$'";
	gpsf+="plot '"+filename_+"-structure-factor.dat' u 1:2 lt 1 lc 6 t 'real',\\";
	gpsf+="     '"+filename_+"-structure-factor.dat' u 1:3 lt 1 lc 7 t 'imag'";
	gpsf.save_file();
	//gpsf.create_image(true);
	rst_file_->link_figure(analyse_+path_+dir_+filename_+"-structure-factor.png","Structure factor",analyse_+path_+dir_+filename_+"-structure-factor.gp",1000);
	/*}*/
	/*!save some additionnal values */
	/*{*/
	jd_write_->write("delta",delta_);
	jd_write_->write("energy per site",E_);
	jd_write_->write("polymerization strength",poly_e(N_/m_-1)-poly_e(N_/m_-2));
	jd_write_->write("critical exponents",exponents);
	/*}*/

	rst_file_->text(read_->get_header());
	rst_file_->save(false);
	delete rst_file_;
	rst_file_ = NULL;

	return tostring(delta_);
}

std::string ChainPolymerized::extract_level_6(){
	Data<double> tmp_E;
	E_.set_x(1e33);
	Vector<double> exponents;
	Vector<double> tmp_exponents;
	double polymerization_strength;
	double tmp_polymerization_strength;
	double tmp_delta;
	unsigned int nof(0);
	(*read_)>>nof;

	bool checked_zero(false);
	for(unsigned int i(0);i<nof;i++){
		(*read_)>>tmp_delta>>tmp_E>>tmp_polymerization_strength>>tmp_exponents;
		if(tmp_E.get_x()<E_.get_x()){ 
			E_ = tmp_E;
			delta_ = tmp_delta;
			polymerization_strength = tmp_polymerization_strength;
			exponents = tmp_exponents;
		}
		if(are_equal(tmp_delta,0)){ checked_zero = true;}
	}

	if(!checked_zero){ std::cerr<<"need to run delta = 0 for "<<N_<<" "<<m_<<" "<<n_<<" "<<bc_<<std::endl; }

	jd_write_->add_to_header("\n");
	save();
	jd_write_->write("energy per site",E_);
	jd_write_->write("polymerization strength",polymerization_strength);
	jd_write_->write("critical exponents",exponents);

	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp+="set title '$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$ bc=$"+tostring(bc_)+"$'";
	gp+="set xlabel '$\\delta$' offset 0,1";
	gp+="set y2label '$\\dfrac{E}{n}$' rotate by 0";
	gp.xrange("0.0","");
	gp+="f(x) = "+std::string(are_equal(delta_,0)?"a+b*x**c":"a+b*(x-c)*(x-c)"); 
	gp+="a="+tostring(E_.get_x());
	gp+="b=1";
	gp+="c=1";
	gp+="set fit quiet";
	gp+="fit f(x) '"+filename_+".dat' u 1:($6==0?$2:1/0):3 zerror via a,b,c";
	gp+="plot '"+filename_+".dat' u 1:(($5==0 && $6==1)?$2:1/0):3 lc 5 w e t 'Not Converged',\\";
	gp+="     '"+filename_+".dat' u 1:(($5==1 && $6==1)?$2:1/0):3 lc 6 w e t 'Converged',\\";
	gp+="     '"+filename_+".dat' u 1:($6==0?$2:1/0):3 lc 7 w e t 'Mean',\\";
	gp+="     f(x) lc 7 lw 0.5 "+std::string(are_equal(delta_,0)?"notitle":"t sprintf('min %3.4f',c)");

	gp.save_file();
	gp.create_image(true);

	return filename_;
}

std::string ChainPolymerized::extract_level_5(){
	double polymerization_strength;
	Vector<double> exponents;
	(*read_)>>E_>>polymerization_strength>>exponents;

	jd_write_->add_to_header("\n");
	save();
	jd_write_->write("energy per site",E_);
	jd_write_->write("polymerization strength",polymerization_strength);
	jd_write_->write("critical exponents",exponents);

	std::cerr<<"level5 "<<N_<<" "<<m_<<" "<<n_<<" "<<bc_<<" "<<delta_<<" "<<E_<<" "<<exponents<<std::endl;;
	return filename_;
}
/*}*/
