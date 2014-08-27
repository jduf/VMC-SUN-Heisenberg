#include "ChainPolymerized.hpp"

ChainPolymerized::ChainPolymerized(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc, double delta):
	System(ref,N,m,n,M,bc),
	Chain<double>(N_/m_,"chain-polymerized"),
	delta_(delta)
{
	if(status_==1){
		init_fermionic();

		filename_ += "-delta" + tostring(delta_);
		system_info_.text("+ Spin chain, with different real hopping term.");
		system_info_.text("  For N colors and m particules per sites, every");
		system_info_.text("  N/m, there is a weaker bound, namely t-delta");
		system_info_.text("  instead of t+delta. (t=1,delta>0) :");
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
	for(unsigned int i(0); i < n_; i += a){
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
	data_write_->precision(10);

	IOFiles corr_file(analyse_+path_+dir_+filename_+"-corr.dat",true);
	IOFiles lr_corr_file(analyse_+path_+dir_+filename_+"-long-range-corr.dat",true);//should not be declared when type!=2

	Vector<double> poly_e(N_/m_,0);
	Vector<double> lrc_mean;
	unsigned int nruns;
	unsigned int tmax;

	(*read_)>>nruns>>tmax;
	(*data_write_)<<"% delta E dE 0|1"<<IOFiles::endl;
	/* the +1 is the averages over all runs */
	for(unsigned int i(0);i<nruns+1;i++){ 
		(*read_)>>E_>>corr_>>lr_corr_;
		if(lrc_mean.size() == 0){ lrc_mean.set(lr_corr_.size(),0);}
		if(i<nruns){
			unsigned int k(0);
			while(k<corr_.size()){
				/*the condition k<corr_.size() is usefull with open boundary conditions*/
				for(unsigned int j(0);j<N_/m_ && k<corr_.size() ;j++){
					poly_e(j) += corr_[k].get_x();
					k++;
				}
			}
		}
		(*data_write_)<<delta_<<" "<<E_.get_x()<<" "<<E_.get_dx()<<" "<<(i<nruns?true:false)<<IOFiles::endl;
		for(unsigned int j(0);j<corr_.size();j++){
			corr_file<<j+0.5<<" "<<corr_[j]<<" "<<(i<nruns?true:false)<<IOFiles::endl;
		}
		for(unsigned int j(0);j<lr_corr_.size();j++){
			lr_corr_file<<j<<" "<<lr_corr_[j]<<" "<<(i<nruns?true:false)<<IOFiles::endl;
			if(i<nruns){ lrc_mean(j) += lr_corr_[j].get_x();}
		}
		if(this->lr_corr_.size()>0){
			lr_corr_file<<n_<<" "<<lr_corr_[0]<<" "<<(i<nruns?true:false)<<IOFiles::endl;
		}
	}
	poly_e /= nruns*n_*m_/N_;
	poly_e.sort(std::less<double>());

	jd_write_->write("delta",delta_);
	jd_write_->write("energy per site",E_);
	jd_write_->write("polymerization strength",poly_e(N_/m_-1)-poly_e(N_/m_-2));

	/*{*/
	Gnuplot gp(analyse_+path_+dir_,filename_+"-corr");
	gp+="set xlabel 'site' offset 0,0.5";
	gp+="set ylabel '$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(i+1)>$' offset 1";
	gp+="set title '$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$ bc="+tostring(bc_)+" $\\delta="+tostring(delta_)+"$'";
	gp+="plot '"+filename_+"-corr.dat' u 1:($6==1?$2:1/0):3 w errorbars lt 1 lc 1 lw 2 t 'Independant measures',\\";
	gp+="     '"+filename_+"-corr.dat' u 1:($6==0?$2:1/0):3 w errorbars lt 1 lc 2 lw 2 t 'Mean',\\";
	gp+="     "+tostring(poly_e(N_/m_-1)) + " w l lc 3 t 'd-merization="+tostring(poly_e(N_/m_-1)-poly_e(N_/m_-2))+"',\\";
	gp+="     "+tostring(poly_e(N_/m_-2)) + " w l lc 3 notitle";
	gp.save_file();
	rst_file_->link_figure(analyse_+path_+dir_+filename_+"-corr.png","Correlation on links",analyse_+path_+dir_+filename_+"-corr.gp",1000);

	if(this->lr_corr_.size()>0){
		unsigned int llr(lr_corr_.size());
		Vector<std::complex<double> > Ck(llr,0.0);
		std::complex<double> normalize(0.0);
		double dk(2.0*M_PI/llr);

		lrc_mean /= nruns;
		for(unsigned int k(0);k<llr;k++){
			for(unsigned int i(0);i<llr;i++){
				Ck(k) += std::polar(lrc_mean(i),dk*k*i);
			}
			normalize += Ck(k); 
		}
		Ck /= dk*normalize;

		IOFiles data_sf(analyse_+path_+dir_+filename_+"-structure-factor.dat",true);
		for(unsigned int k(0);k<llr;k++){
			data_sf<<dk*k<<" "<<Ck(k).real()<<" "<<Ck(k).imag()<<" "<<std::abs(Ck(k))<<IOFiles::endl;
		}

		Gnuplot gplr(analyse_+path_+dir_,filename_+"-long-range-corr");
		gplr.xrange(0,llr+1);
		gplr+="set xlabel '$\\|i-j\\|$' offset 0,0.5";
		gplr+="set ylabel '$<S_{\\alpha}^{\\alpha}(i)S_{\\alpha}^{\\alpha}(j)>-\\dfrac{m}{N}$' offset 1";
		gplr+="set title '$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$ bc="+tostring(bc_)+" $\\delta="+tostring(delta_)+"$'";
		gplr+="set key right bottom";
		gplr+="m="+tostring(m_)+".0";
		gplr+="N="+tostring(N_)+".0";
		gplr+="n="+tostring(n_)+".0";
		gplr+="f(x) = a/(x*x) + b*cos(2.0*pi*x*m/N)/(x**eta)";
		gplr+="f(x) = bpi*cos(2.0*pi*x*m/N)*(x**(-ap)+(n-x)**(-ap))+b0*(x**(-a0)+(n-x)**(-a0))";
		gplr+="set fit quiet";
		gplr+="fit ["+tostring(spuc_*ceil(L_*0.2))+":"+tostring(spuc_*floor(L_*0.8))+"] f(x) '"+filename_+"-long-range-corr.dat' u 1:($6==0?$2:1/0) noerrors via b0,bpi,ap,a0"; 
		gplr+="plot '"+filename_+"-long-range-corr.dat' u 1:($6==1?$2:1/0):3 w errorbars lt 1 lc 1 lw 2 t 'Independant measures',\\";
		gplr+="     '"+filename_+"-long-range-corr.dat' u 1:($6==0?$2:1/0):3 w errorbars lt 1 lc 2 lw 2 t 'Mean',\\";
		gplr+="     f(x) t sprintf('$a_\\pi=%f$, $a_0=%f$',ap,a0)";
		gplr.save_file();
		gplr.create_image(true);
		rst_file_->link_figure(analyse_+path_+dir_+filename_+"-long-range-corr.png","Long range correlation",analyse_+path_+dir_+filename_+"-long-range-corr.gp",1000);

		Gnuplot gpsf(analyse_+path_+dir_,filename_+"-structure-factor");
		gpsf+="set title '$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$ bc="+tostring(bc_)+" $\\delta="+tostring(delta_)+"$'";
		gpsf+="set key bottom";
		gpsf.xrange("0","2*pi");
		switch(N_/m_){
			case 3: { gpsf+="set xtics ('0' 0,'$2\\pi/3$' 2.0*pi/3.0,'$\\pi$' pi,'$4\\pi/3$' 4.0*pi/3.0,'$2\\pi$' 2.0*pi)"; } break;
			case 5: { gpsf+="set xtics ('0' 0,'$2\\pi/5$' 2.0*pi/5.0, '$4\\pi/5$' 4.0*pi/5.0, '$6\\pi/5$' 6.0*pi/5.0, '$8\\pi/5$' 8.0*pi/5.0, '$2\\pi$' 2.0*pi)"; } break;
			default:{ gpsf+="set xtics ('0' 0,'$\\pi/2$' pi/2.0,'$\\pi$' pi,'$3\\pi/2$' 3.0*pi/2.0,'$2\\pi$' 2.0*pi)"; } break;
		}
		gpsf+="plot '"+filename_+"-structure-factor.dat' u 1:2 t 'real',\\";
		gpsf+="     '"+filename_+"-structure-factor.dat' u 1:3 t 'imag',\\";
		gpsf+="     '"+filename_+"-structure-factor.dat' u 1:4 t 'norm'";
		gpsf.save_file();
		gpsf.create_image(true);
		rst_file_->link_figure(analyse_+path_+dir_+filename_+"-structure-factor.png","Structure factor",analyse_+path_+dir_+filename_+"-structure-factor.gp",1000);
	}
	/*}*/

	rst_file_->text(read_->get_header());
	rst_file_->save(false);
	delete rst_file_;
	rst_file_ = NULL;

	return tostring(delta_);
}

std::string ChainPolymerized::extract_level_6(){
	double min_delta(delta_);
	double min_polymerization_strength(0.0);
	double polymerization_strength;
	Data<double> min_E;
	min_E.set_x(0.0);

	unsigned int nof(0);
	unsigned int idx(0);
	(*read_)>>nof;
	for(unsigned int i(0);i<nof;i++){
		(*read_)>>delta_>>E_>>polymerization_strength;
		if(E_.get_x()<min_E.get_x()){ 
			idx = i;
			min_E = E_;
			min_delta = delta_;
			min_polymerization_strength = polymerization_strength;
		}
	}
	delta_ = min_delta;

	jd_write_->add_to_header("\n");
	save();
	jd_write_->write("energy per site",min_E);
	jd_write_->write("polymerization strength",min_polymerization_strength);

	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp+="set xlabel '$\\delta$' offset 0,1";
	gp+="set ylabel '$\\dfrac{E}{n}$' rotate by 0 offset 1";
	if(idx==0){
		gp+="f(x) = a+b*x**eta";
		gp+="a="+tostring(min_E.get_x());
		gp+="b=1";
		gp+="eta=1";
		gp+="set fit quiet";
		gp+="fit f(x) '"+filename_+".dat' u 1:($4==0?$2:1/0):3 zerror via a,b,eta";
		gp+="set title '$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$'";
		gp+="plot '"+filename_+".dat' u 1:($4==1?$2:1/0):3 w e t 'Independant measures',\\";
		gp+="     '"+filename_+".dat' u 1:($4==0?$2:1/0):3 w e t 'Mean',\\";
		gp+="     f(x) t sprintf('eta %3.4f',eta)";
	} else {
		gp+="f(x) = a+b*(x-c)*(x-c)";
		gp+="a="+tostring(min_E.get_x());
		gp+="b=1";
		gp+="c="+tostring(delta_);
		gp+="set fit quiet";
		gp+="fit f(x) '"+filename_+".dat' u 1:($4==0?$2:1/0):3 zerror via a,b,c";
		gp+="set title '$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$'";
		gp+="plot '"+filename_+".dat' u 1:($4==1?$2:1/0):3 w e t 'Independant measures',\\";
		gp+="     '"+filename_+".dat' u 1:($4==0?$2:1/0):3 w e t 'Mean',\\";
		gp+="     f(x) t sprintf('min %3.4f',c)";
	}
	gp.save_file();
	gp.create_image(true);

	return filename_;
}
/*}*/
