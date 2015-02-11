#include "ChainPolymerized.hpp"

ChainPolymerized::ChainPolymerized(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc, Vector<double> ti):
	System(ref,N,m,n,M,bc),
	Chain<double>((are_equal(ti,Vector<double>(N_/m_,1.0))?1:N_/m_),"chain-polymerized"),
	ti_(ti)
{
	if(status_==2){
		init_fermionic();
		filename_ += "-ti";
		for(unsigned int j(0);j<ti_.size();j++){
			filename_ += ((ti_(j)>0)?"+":"")+tostring(ti_(j));
		}
		if(spuc_ != 1){
			system_info_.text("+ Spin chain, with different real hopping terms.");
			if(spuc_ != 4){
				std::string tmp("");
				for(unsigned int i(0);i<spuc_-1;i++){ tmp += "=\\bullet"; }
				system_info_.text("+  :math:`t_i : "+tmp+"-`");
			}
			else{system_info_.text("+  :math:`t_i : \\equiv\\bullet=\\bullet\\equiv\\bullet-`");}
		} else {
			system_info_.text("+ Spin chain with real and uniform hopping terms");
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
			H_(i+j,nb(0,0)) = nb(0,1)*ti_(j);
		}
	}
	H_ += H_.transpose();
}

void ChainPolymerized::create(){
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);
	lr_corr_.set(links_.row(),50,5,false);

	compute_H();
	diagonalize(true);
	for(unsigned int c(0);c<N_;c++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_(c);j++){
				EVec_[c](i,j) = H_(i,j);
			}
		}
	}
	if(status_==2){
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
			n1 += norm_squared(tmp1);
			n2 += norm_squared(tmp2);
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

void ChainPolymerized::save() const {
	GenericSystem<double>::save();
	if(spuc_==4){jd_write_->write("ti ("+tostring(ti_(1))+","+tostring(ti_(3))+")",ti_);}
	else{jd_write_->write("ti ("+tostring(ti_(spuc_-1))+")",ti_);}
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
	std::string ti_string("(");
	for(unsigned int j(0);j<ti_.size()-1;j++){
		ti_string += tostring(ti_(j))+",";
	}
	ti_string += tostring(ti_(ti_.size()-1))+")";
	std::string title("$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$ bc="+tostring(bc_)+" $t_{ij}="+ti_string+"$");

	/*!extract jdbin*/
	/*{*/
	IOFiles corr_file(analyse_+path_+dir_+filename_+"-corr.dat",true);
	IOFiles lr_corr_file(analyse_+path_+dir_+filename_+"-long-range-corr.dat",true);

	Vector<double> lrc_mean(links_.row(),0);
	Vector<double> poly_e(N_/m_,0);
	unsigned int nruns;
	unsigned int tmax;

	(*read_)>>nruns>>tmax;
	(*data_write_)<<"% ti E dx conv(0|1) #conv mean(0|1)"<<IOFiles::endl;
	corr_file<<"%(2i+1)/2 corr(i,i+1) dx conv(0|1) #conv mean(0|1)"<<IOFiles::endl;
	lr_corr_file<<"%j corr(i,j) dx conv(0|1) #conv mean(0|1)"<<IOFiles::endl;
	/*!the +1 is the average over all runs */
	for(unsigned int i(0);i<nruns+1;i++){ 
		(*read_)>>E_>>corr_>>lr_corr_;
		(*data_write_)<<ti_<<" "<<E_<<" "<<(i<nruns)<<IOFiles::endl;
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
	gp+="set key center";
	gp.label("x","site","offset 0,0.5");
	gp.label("y2","$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(i+1)>$");
	gp.title(title);
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
	unsigned int xi;
	unsigned int xf;
	Vector<double> exponents;
	bool fit(compute_critical_exponents(lrc_mean,xi,xf,exponents));

	Gnuplot gplr(analyse_+path_+dir_,filename_+"-long-range-corr");
	gplr.range("x",N_/m_,n_-N_/m_);
	gplr.label("x","$\\|i-j\\|$","offset 0,0.5");
	gplr.label("y2","$<S_{\\alpha}^{\\alpha}(i)S_{\\alpha}^{\\alpha}(j)>-\\dfrac{m^2}{N}$","offset 1");
	gp.title(title);
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
	gplr+="     f(x) lc 7 " + std::string(fit?"lw 0.5":"dt 2") + " t sprintf('$\\eta=%f$, $\\mu=%f$',p1,p3)";
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
	rst_file_->link_figure(analyse_+path_+dir_+filename_+"-structure-factor.png","Structure factor",analyse_+path_+dir_+filename_+"-structure-factor.gp",1000);
	/*}*/
	/*!save some additionnal values */
	/*{*/
	if(spuc_==4){jd_write_->write("ti ("+tostring(ti_(1))+","+tostring(ti_(3))+")",ti_);}
	else{jd_write_->write("ti ("+tostring(ti_(spuc_-1))+")",ti_);}
	jd_write_->write("energy per site",E_);
	jd_write_->write("polymerization strength",poly_e(N_/m_-1)-poly_e(N_/m_-2));
	jd_write_->write("critical exponents",exponents);
	/*}*/

	rst_file_->text(read_->get_header());
	rst_file_->save(false);
	delete rst_file_;
	rst_file_ = NULL;

	return ti_string;
}

std::string ChainPolymerized::extract_level_6(){
	Data<double> tmp_E;
	E_.set_x(1e33);
	Vector<double> exponents;
	Vector<double> tmp_exponents;
	double polymerization_strength;
	double tmp_polymerization_strength;
	Vector<double> tmp_ti;
	unsigned int nof(0);
	(*read_)>>nof;

	for(unsigned int i(0);i<nof;i++){
		(*read_)>>tmp_ti>>tmp_E>>tmp_polymerization_strength>>tmp_exponents;
		if(tmp_E.get_x()<E_.get_x()){ 
			E_ = tmp_E;
			ti_ = tmp_ti;
			polymerization_strength = tmp_polymerization_strength;
			exponents = tmp_exponents;
		}
	}

	jd_write_->add_header()->nl();
	save();
	jd_write_->write("energy per site",E_);
	jd_write_->write("polymerization strength",polymerization_strength);
	jd_write_->write("critical exponents",exponents);

	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp.title("$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$ bc=$"+tostring(bc_)+"$");
	if(N_/m_!=4){
		gp.label("x","$t_"+tostring(N_/m_)+"$"," offset 0,1");
		gp.label("y2","$\\dfrac{E}{n}$","rotate by 0");
		gp.range("x","0.0","");
		gp+="f(x) = "+std::string(are_equal(ti_(N_/m_-1),0)?"a+b*x**c":"a+b*(x-c)*(x-c)"); 
		gp+="a="+tostring(E_.get_x());
		gp+="b=1";
		gp+="c=1";
		gp+="set fit quiet";
		gp+="fit f(x) '"+filename_+".dat' u "+tostring(N_/m_)+":($" +tostring(N_/m_+5)+"==0?$"                              +tostring(N_/m_+1)+":1/0):"+tostring(N_/m_+2)+" zerror via a,b,c";
		gp+="plot '"+filename_+".dat' u "    +tostring(N_/m_)+":(($"+tostring(N_/m_+4)+"==0 && $"+tostring(N_/m_+5)+"==1)?$"+tostring(N_/m_+1)+":1/0):"+tostring(N_/m_+2)+" lc 5 w e t 'Not Converged',\\";
		gp+="     '"+filename_+".dat' u "    +tostring(N_/m_)+":(($"+tostring(N_/m_+4)+"==1 && $"+tostring(N_/m_+5)+"==1)?$"+tostring(N_/m_+1)+":1/0):"+tostring(N_/m_+2)+" lc 6 w e t 'Converged',\\";
		gp+="     '"+filename_+".dat' u "    +tostring(N_/m_)+":($" +tostring(N_/m_+5)+"==0?$"                              +tostring(N_/m_+1)+":1/0):"+tostring(N_/m_+2)+" lc 7 w e t 'Mean',\\";
		gp+="     f(x) lc 7 lw 0.5 "+std::string(are_equal(ti_(N_/m_-1),0)?"notitle":"t sprintf('min %3.4f',c)");
	} else {
		gp.label("x","$t_2$");
		gp.label("y","$t_4$");
		gp.label("z","$\\dfrac{E}{n}$","rotate by 0");
		gp+="set xyplane 0";
		gp+="splot '"+filename_+".dat' u 2:4:(($8==0 && $9==1)?$5:1/0) lc 5 t 'Not Converged',\\";
		gp+="      '"+filename_+".dat' u 2:4:(($8==1 && $9==1)?$5:1/0) lc 6 t 'Converged',\\";
		gp+="      '"+filename_+".dat' u 2:4:($9==0?$5:1/0) lc 7 t 'Mean'";
	}

	gp.save_file();
	gp.create_image(true);

	return filename_;
}
/*}*/
