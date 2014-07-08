#include"ChainFermi.hpp"

ChainFermi::ChainFermi(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc):
	System(ref,N,m,n,M,bc),
	Chain<double>(1,"chain-fermi")
{
	if(status_==1){
		init_fermionic();
		compute_T();

		system_info_.text("Spin ChainFermi, all the hopping parameters are real");
	}
}

/*{method needed for running*/
void ChainFermi::compute_T(){
	double t(1.0);
	T_.set(n_,n_,0);
	Matrix<int> nb;
	for(unsigned int i(0); i< n_; i++){
		nb = get_neighbourg(i);
		T_(i,nb(0,0)) = nb(0,1)*t;
	}
	T_ += T_.transpose();
}

void ChainFermi::create(){
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);
	//if(type==2){ long_range_corr_.set(n_/3); }

	diagonalize_T();
	for(unsigned int c(0);c<N_;c++){
		if(!is_degenerate(c)){
			EVec_[c].set(n_,M_(c));
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = T_(i,j);
				}
			}
		}
	}
}
/*}*/

/*{method needed for checking*/
void ChainFermi::check(){
	BandStructure<double> bs(T_,Lx_,spuc_,bc_);
}
/*}*/

/*{method needed for analysing*/
std::string ChainFermi::extract_level_6(){
	rst_file_ = new RSTFile(info_+path_+dir_,filename_);

	unsigned int nruns;
	unsigned int tmax;

	(*read_)>>nruns>>tmax;
	//std::cout<<delta_<<" "<<N_<<" "<<m_<<" "<<n_<<" "<<M_<<" "<<tmax<<" "<<nruns<<std::endl;
	IOFiles corr_file(analyse_+path_+dir_+filename_+"-corr.dat",true);
	IOFiles long_range_corr_file(analyse_+path_+dir_+filename_+"-long-range-corr.dat",true);//should not be delcared when type!=2
	data_write_->precision(10);
	(*data_write_)<<"% E dE 0|1"<<IOFiles::endl;
	/* the +1 is the averages over all runs */
	Vector<double> poly_e(N_/m_,0);
	for(unsigned int i(0);i<nruns+1;i++){ 
		(*read_)>>E_>>corr_>>long_range_corr_;
		if(i<nruns){
			unsigned int k(0);
			while(k<corr_.size()){
				for(unsigned int j(0);j<N_/m_;j++){
					poly_e(j) += corr_[k].get_x();
					k++;
				}
			}
		}
		(*data_write_)<<" "<<E_.get_x()<<" "<<E_.get_dx()<<" "<<(i<nruns?true:false)<<IOFiles::endl;
		for(unsigned int j(0);j<corr_.size();j++){
			corr_file<<j+0.5<<" "<<corr_[j]<<" "<<(i<nruns?true:false)<<IOFiles::endl;
		}
		for(unsigned int j(0);j<long_range_corr_.size();j++){
			long_range_corr_file<<j+1<<" "<<long_range_corr_[j]<<" "<<(i<nruns?true:false)<<IOFiles::endl;
		}

	}
	poly_e /= nruns*n_*m_/N_;
	poly_e.sort(std::less<double>());

	(*jd_write_)("E",E_);
	(*jd_write_)("polymerization strength",poly_e(N_/m_-1)-poly_e(N_/m_-2));

	/*{*/
	Gnuplot gp(analyse_+path_+dir_,filename_+"-corr");
	gp+="set xlabel 'site' offset 0,0.5";
	gp+="set ylabel '$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(i+1)>$' offset 1";
	gp+="set title '$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$ bc="+tostring(bc_)+"'";
	gp+="plot '"+filename_+"-corr.dat' u 1:($6==1?$2:1/0):3 w errorbars lt 1 lc 1 lw 2 t 'Independant measures',\\";
	gp+="     '"+filename_+"-corr.dat' u 1:($6==0?$2:1/0):3 w errorbars lt 1 lc 2 lw 2 t 'Mean',\\";
	gp+="     "+tostring(poly_e(N_/m_-1)) + " w l lc 3 t 'd-merization="+tostring(poly_e(N_/m_-1)-poly_e(N_/m_-2))+"',\\";
	gp+="     "+tostring(poly_e(N_/m_-2)) + " w l lc 3 notitle";
	gp.save_file();
	rst_file_->link_figure(analyse_+path_+dir_+filename_+"-corr.png","Correlation on links",analyse_+path_+dir_+filename_+"-corr.gp",1000);

	unsigned int length(long_range_corr_.size());
	if(length>0){
		Gnuplot gp(analyse_+path_+dir_,filename_+"-long-range-corr");
		gp+="stats '"+filename_+"-long-range-corr.dat' nooutput";
		gp.xrange(0,length+1);
		gp.yrange("1.1*STATS_min_y","1.1*STATS_max_y");
		gp+="set xlabel '$\\|i-j\\|$' offset 0,0.5";
		gp+="set ylabel '$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(j)>$' offset 1";
		gp+="set title '$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$ bc="+tostring(bc_)+"'";
		gp+="set key right bottom";
		gp+="a=1.0";
		gp+="b=1.0";
		gp+="eta=1.0";
		gp+="m="+tostring(m_)+".0";
		gp+="N="+tostring(N_)+".0";
		gp+="f(x) = a/(x*x) + b*cos(2.0*pi*x*m/N)/(x**eta)";
		gp+="set fit quiet";
		switch(N_/m_){
			case 2:{ gp+="fit [3:"+tostring(length)+"] f(x) '"+filename_+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; } break;
			case 3:{
					   switch((length + 1) % 3){
						   case 0:{ gp+="fit [2:"+tostring(length)+"] f(x) '"+filename_+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
						   case 1:{ gp+="fit [5:"+tostring(length)+"] f(x) '"+filename_+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
						   case 2:{ gp+="fit [3:"+tostring(length)+"] f(x) '"+filename_+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
					   }break;
				   }break;
			default :{ gp+="fit ["+tostring(N_-1)+":] f(x) '"+filename_+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
		}
		gp+="plot for [IDX=0:"+tostring(nruns-1)+"] '"+filename_+"-long-range-corr.dat' i IDX u 1:($4==1?$2:1/0):3 w errorbars lt 1 lc 3 ps 0 notitle,\\";
		gp+="     for [IDX=0:"+tostring(nruns-1)+"] '"+filename_+"-long-range-corr.dat' i IDX u 1:($4==0?$2:1/0):3 w errorbars lt 1 lc 1 ps 0 notitle,\\";
		gp+="                   '"+filename_+"-long-range-corr.dat' i "+tostring(nruns)+" u 1:2:3 w errorbars lt 1 lc 2 lw 2 notitle,\\";
		gp+="                   f(x) notitle";
		gp.save_file();
		rst_file_->link_figure(analyse_+path_+dir_+filename_+"-long-range-corr.png","Long range correlation",analyse_+path_+dir_+filename_+"-long-range-corr.gp",1000);
	}
	/*}*/

	rst_file_->text(read_->get_header());
	rst_file_->save(false);
	delete rst_file_;
	rst_file_ = NULL;

	return filename_;
}

std::string ChainFermi::extract_level_5(){
	double min_polymerization_strength(0.0);
	double polymerization_strength;
	Data<double> min_E;
	min_E.set_x(0.0);

	unsigned int nof(0);
	(*read_)>>nof;
	for(unsigned int i(0);i<nof;i++){
		(*read_)>>E_>>polymerization_strength;
		if(E_.get_x()<min_E.get_x()){ 
			min_E = E_;
			min_polymerization_strength = polymerization_strength;
		}
	}

	save(*jd_write_);
	(*jd_write_)("energy per site",min_E);
	(*jd_write_)("polymerization strength",min_polymerization_strength);

	return filename_;
}

std::string ChainFermi::extract_level_3(){
	double polymerization_strength;
	(*read_)>>E_>>polymerization_strength;
	(*data_write_)<<n_<<" "<<polymerization_strength<<IOFiles::endl;

	return filename_;
}
/*}*/
