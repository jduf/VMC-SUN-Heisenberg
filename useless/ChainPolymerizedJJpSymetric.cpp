#include "ChainPolymerizedJJp.hpp"

ChainPolymerizedJJp::ChainPolymerizedJJp(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc, Vector<double> const& t, Vector<double> const& JJp):
	System(ref,N,m,n,M,bc),
	Chain<double>((my::are_equal(t,Vector<double>(N_/m_,1.0))?1:N_/m_),"chain-polymerized-jjp"),
	t_(t),
	JJp_(JJp)
{
	if(status_==2){
		init_fermionic();
		filename_ += "-t";
		for(unsigned int j(0);j<t_.size();j++){
			filename_ += ((t_(j)>0)?"+":"")+my::tostring(t_(j));
		}
		filename_ += "-J";
		for(unsigned int j(0);j<JJp_.size();j++){
			filename_ += ((JJp_(j)>0)?"+":"")+my::tostring(JJp_(j));
		}
		for(unsigned int i(0);i<J_.size();i++){
			if(i%2){ J_(i) = JJp(1); }
			else { J_(i) = JJp_(0); }
		}
		system_info_.text("The Hamiltonian is of the form : ");
		system_info_.math_centered("H=JS_iS_{i+1} + J'S_{i+1}S_{i+2}");
		system_info_.text("with " + RST::math("J<0") + " and " +RST::math("J'<0"));
	}
}

/*{method needed for running*/
void ChainPolymerizedJJp::compute_H(unsigned int const& c){
	H_.set(n_,n_,0);
	Matrix<int> nb;
	//int bc(1);
	//unsigned int s;
	//for(unsigned int i(0);i<n_;i++){
		//nb = get_neighbourg(i);
		//bc = nb(0,1);
		//if(i%2==0){ H_(i,nb(0,0)) = (c==1?-1.0:1.0)*bc*t_(0); }
		//else { H_(i,nb(0,0)) = bc*t_(1);  }
		//s = nb(0,0);
		//nb = get_neighbourg(s);
		//bc *= nb(0,1);
		//H_(i,nb(0,0)) = bc*t_(1);
		//s = nb(0,0);
		//nb = get_neighbourg(s);
		//bc *= nb(0,1);
		//H_(i,nb(0,0)) = (c==1?-1.0:1.0)*bc*t_(1);
	//}
	
	//for(unsigned int i(0);i<n_;i+=2){
		//nb = get_neighbourg(i);
		//H_(i,nb(0,0)) = (c==1?-1.0:1.0)*nb(0,1)*t_(0);
		//nb = get_neighbourg(i+1);
		//H_(i+1,nb(0,0)) = nb(0,1)*t_(1);
	//}
	
	if(c){
		for(unsigned int i(0);i<n_;i+=2){
			nb = get_neighbourg(i);
			H_(i,nb(0,0)) = -nb(0,1)*t_(0);
			nb = get_neighbourg(i+1);
			H_(i+1,nb(0,0)) = nb(0,1)*t_(1);
		}
	} else {
		for(unsigned int i(0);i<n_;i+=2){
			nb = get_neighbourg(i);
			H_(i,nb(0,0)) = nb(0,1)*t_(0);
			nb = get_neighbourg(i+1);
			H_(i+1,nb(0,0)) = nb(0,1)*t_(1);
		}
	}

	H_ += H_.transpose();
}

void ChainPolymerizedJJp::create(){
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);
	lr_corr_.set(links_.row(),50,5,false);

	for(unsigned int c(0);c<N_;c++){
		compute_H(c);
		diagonalize(true);
		if(status_==1){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = H_(i,j);
				}
			}
			std::cout<<"U"<<std::endl;
			std::cout<<H_<<std::endl;
			if(c<N_-1){ status_++; }
		}
	}
}

void ChainPolymerizedJJp::save() const {
	GenericSystem<double>::save();
	std::string s("");
	for(unsigned int i(0);i<t_.size()-1;i++){ s += my::tostring(t_(i))+","; }
	s += my::tostring(t_(t_.size()-1));
	jd_write_->write("t ("+s+")",t_);

	s="";
	for(unsigned int i(0);i<JJp_.size()-1;i++){ s += my::tostring(JJp_(i))+","; }
	s += my::tostring(JJp_(JJp_.size()-1));
	jd_write_->write("JJp ("+s+")",JJp_);
}
/*}*/

/*{method needed for checking*/
void ChainPolymerizedJJp::check(){
	for(unsigned int i(0);i<J_.size();i++){
		std::cout<<J_(i)<<" "<<links_(i,0)<<" "<<links_(i,1)<<std::endl;
	}
	compute_H(0);
	std::cout<<"H(0)"<<std::endl;
	std::cout<<H_<<std::endl;
	compute_H(1);
	std::cout<<"H(1)"<<std::endl;
	std::cout<<H_<<std::endl;
	create();
	std::cout<<"U(0)"<<std::endl;
	std::cout<<EVec_[0]<<std::endl;
	std::cout<<"U(1)"<<std::endl;
	std::cout<<EVec_[1]<<std::endl;
}
/*}*/

/*{method needed for analysing*/
std::string ChainPolymerizedJJp::extract_level_7(){
	rst_file_ = new RSTFile(info_+path_+dir_,filename_);
	std::string t_string("(");
	for(unsigned int i(0);i<JJp_.size()-1;i++){
		t_string += my::tostring(JJp_(i))+",";
	}
	t_string += my::tostring(JJp_(JJp_.size()-1))+"|";
	for(unsigned int i(0);i<t_.size()-1;i++){
		t_string += my::tostring(t_(i))+",";
	}
	t_string += my::tostring(t_(t_.size()-1))+")";
	std::string title("$N="+my::tostring(N_)+"$ $m="+my::tostring(m_)+"$ $n="+my::tostring(n_)+"$ bc="+my::tostring(bc_)+" $t_{ij}="+t_string+"$");

	/*!extract jdbin*/
	/*{*/
	IOFiles corr_file(analyse_+path_+dir_+filename_+"-corr.dat",true);
	IOFiles lr_corr_file(analyse_+path_+dir_+filename_+"-long-range-corr.dat",true);

	Vector<double> lrc_mean(links_.row(),0);
	Vector<double> poly_e(N_/m_,0);
	unsigned int nruns;
	unsigned int tmax;

	(*read_)>>nruns>>tmax;
	(*data_write_)<<"% J t E dx conv(0|1) #conv mean(0|1)"<<IOFiles::endl;
	corr_file<<"%(2i+1)/2 corr(i,i+1) dx conv(0|1) #conv mean(0|1)"<<IOFiles::endl;
	lr_corr_file<<"%j corr(i,j) dx conv(0|1) #conv mean(0|1)"<<IOFiles::endl;
	/*!the +1 is the average over all runs */
	for(unsigned int i(0);i<nruns+1;i++){
		(*read_)>>E_>>corr_>>lr_corr_;
		(*data_write_)<<JJp_<<t_<<E_<<" "<<(i<nruns)<<IOFiles::endl;
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
	/*!nearest neighbourg correlatons*/
	/*{*/
	Gnuplot gp(analyse_+path_+dir_,filename_+"-corr");
	gp+="set key center";
	gp.label("x","site","offset 0,0.5");
	gp.label("y2","$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(i+1)>$");
	gp.title(title);
	gp+="plot '"+filename_+"-corr.dat' u 1:(($6==1 && $5==0)?$2:1/0):3 w errorbars lt 1 lc 5 t 'Not converged',\\";
	gp+="     '"+filename_+"-corr.dat' u 1:(($6==1 && $5==1)?$2:1/0):3 w errorbars lt 1 lc 6 t 'Converged',\\";
	gp+="     '"+filename_+"-corr.dat' u 1:($6==0?$2:1/0):3 w errorbars lt 1 lc 7 t 'Mean',\\";
	gp+="     "+my::tostring(poly_e(N_/m_-1)) + " w l lc 3 t 'd-merizaton="+my::tostring(poly_e(N_/m_-1)-poly_e(N_/m_-2))+"',\\";
	gp+="     "+my::tostring(poly_e(N_/m_-2)) + " w l lc 3 notitle";
	gp.save_file();
	//gp.create_image(true);
	rst_file_->link_figure(analyse_+path_+dir_+filename_+"-corr.png","Correlaton on links",analyse_+path_+dir_+filename_+"-corr.gp",1000);
	/*}*/
	/*!long range correlatons*/
	/*{*/
	unsigned int xi;
	unsigned int xf;
	Vector<double> exponents;
	//bool fit(compute_critical_exponents(lrc_mean,xi,xf,exponents));

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
	gplr+="fit [" + my::tostring(xi) + ":" + my::tostring(xf) + "] f(x) '"+filename_+"-long-range-corr.dat' u 1:($6==0?$2:1/0) noerrors via p0,p1,p2,p3";
	gplr+="plot '"+filename_+"-long-range-corr.dat' u 1:(($6==1 && $5==0)?$2:1/0):3 w errorbars lt 1 lc 5 t 'Not converged',\\";
	gplr+="     '"+filename_+"-long-range-corr.dat' u 1:(($6==1 && $5==1)?$2:1/0):3 w errorbars lt 1 lc 6 t 'Converged',\\";
	gplr+="     '"+filename_+"-long-range-corr.dat' u 1:($6==0?$2:1/0):3 w errorbars lt 1 lc 7 t 'Mean'";
	//gplr+="     f(x) lc 7 " + std::string(fit?"lw 0.5":"dt 2") + " t sprintf('$\\eta=%f$, $\\mu=%f$',p1,p3)";
	gplr.save_file();
	//gplr.create_image(true);
	rst_file_->link_figure(analyse_+path_+dir_+filename_+"-long-range-corr.png","Long range correlaton",analyse_+path_+dir_+filename_+"-long-range-corr.gp",1000);
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
		case 3: { gpsf+="set xitcs ('0' 0,'$2\\pi/3$' 2.0*pi/3.0, '$4\\pi/3$' 4.0*pi/3.0,'$2\\pi$' 2.0*pi)"; } break;
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
	/*!save some additonnal values */
	/*{*/
	if(spuc_==4){jd_write_->write("t ("+my::tostring(t_(1))+","+my::tostring(t_(3))+")",t_);}
	else{jd_write_->write("t ("+my::tostring(t_(spuc_-1))+")",t_);}
	jd_write_->write("JJp ("+my::tostring(JJp_(0))+","+my::tostring(JJp_(1))+")",JJp_);
	jd_write_->write("energy per site",E_);
	jd_write_->write("polymerizaton strength",poly_e(N_/m_-1)-poly_e(N_/m_-2));
	jd_write_->write("critical exponents",exponents);
	/*}*/

	rst_file_->text(read_->get_header());
	rst_file_->save(false);
	delete rst_file_;
	rst_file_ = NULL;

	return t_string;
}

std::string ChainPolymerizedJJp::extract_level_6(){
	Data<double> tmp_E;
	E_.set_x(1e33);
	Vector<double> exponents;
	Vector<double> tmp_exponents;
	double polymerizaton_strength;
	double tmp_polymerizaton_strength;
	Vector<double> tmp_t;
	Vector<double> tmp_JJp;
	unsigned int nof(0);
	(*read_)>>nof;

	for(unsigned int i(0);i<nof;i++){
		(*read_)>>tmp_t>>tmp_JJp>>tmp_E>>tmp_polymerizaton_strength>>tmp_exponents;
		if(tmp_E.get_x()<E_.get_x()){
			E_ = tmp_E;
			t_ = tmp_t;
			JJp_ = tmp_JJp;
			polymerizaton_strength = tmp_polymerizaton_strength;
			exponents = tmp_exponents;
		}
	}

	save();
	jd_write_->write("energy per site",E_);
	jd_write_->write("polymerizaton strength",polymerizaton_strength);
	jd_write_->write("critical exponents",exponents);

	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp.title("$N="+my::tostring(N_)+"$ $m="+my::tostring(m_)+"$ $n="+my::tostring(n_)+"$ bc=$"+my::tostring(bc_)+"$ $J=-1$");
	gp.label("x","$J^{\\prime}/|J|$"," offset 0,1");
	gp.label("y2","$\\dfrac{E}{n}$","rotate by 0");
	gp.range("x","0.0","");
	gp+="f(x) = "+std::string(my::are_equal(t_(N_/m_-1),0)?"a+b*x**c":"a+b*(x-c)*(x-c)");
	gp+="a="+my::tostring(E_.get_x());
	gp+="b=1";
	gp+="c=1";
	gp+="set fit quiet";
	std::string c1(my::tostring(JJp_.size()));
	std::string c2(my::tostring(JJp_.size()+t_.size()));
	std::string c3(my::tostring(JJp_.size()+t_.size()+1));//E
	std::string c4(my::tostring(JJp_.size()+t_.size()+2));//dE
	std::string c5(my::tostring(JJp_.size()+t_.size()+4));//conv
	std::string c6(my::tostring(JJp_.size()+t_.size()+5));//mean
	gp+="fit f(x) '"+filename_+".dat' u "+c1+":(($"+c6+"==0)?         $"     +c3+":1/0):"+c4+" zerror via a,b,c";
	gp+="plot     '"+filename_+".dat' u "+c1+":(($"+c5+"==0 && $"+c6+"==1)?$"+c3+":1/0):"+c4+" lc 5 w e t 'Not Converged',\\";
	gp+="         '"+filename_+".dat' u "+c1+":(($"+c5+"==1 && $"+c6+"==1)?$"+c3+":1/0):"+c4+" lc 6 w e t 'Converged',\\";
	gp+="         '"+filename_+".dat' u "+c1+":(($"+c6+"==0)?         $"     +c3+":1/0):"+c4+" lc 7 w e t 'Mean',\\";
	gp+="         f(x) lc 7 lw 0.5 "+std::string(my::are_equal(t_(N_/m_-1),0)?"notitle":"t sprintf('min %3.4f',c)");

	gp.save_file();
	gp.create_image(true);

	return filename_;
}
/*}*/
