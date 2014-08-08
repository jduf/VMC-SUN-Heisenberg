#ifndef DEF_CHAINFERMI
#define DEF_CHAINFERMI

#include "Chain.hpp"

template<typename Type>
class ChainFermi: public Chain<Type>{
	public:
		ChainFermi(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc);
		~ChainFermi(){}

		void create();
		void check();

	private:
		void compute_H();
		std::string extract_level_7();
		std::string extract_level_6();
		std::string extract_level_3();
};

template<typename Type>
ChainFermi<Type>::ChainFermi(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc):
	System(ref,N,m,n,M,bc),
	Chain<Type>(1,"chain-fermi")
{
	if(this->status_==1){
		this->init_fermionic();

		this->system_info_.text("Spin ChainFermi, all the hopping parameters are real");
	}
}

/*{method needed for running*/
template<typename Type>
void ChainFermi<Type>::compute_H(){
	double t(1.0);
	this->H_.set(this->n_,this->n_,0);
	Matrix<int> nb;
	for(unsigned int i(0);i<this->n_;i++){
		nb = this->get_neighbourg(i);
		this->H_(i,nb(0,0)) = nb(0,1)*t;
	}
	this->H_ += this->H_.transpose();
}
/*}*/

/*{method needed for checking*/
template<typename Type>
void ChainFermi<Type>::check(){
	this->create();
}
/*}*/

/*{method needed for analysing*/
template<typename Type>
std::string ChainFermi<Type>::extract_level_7(){
	this->rst_file_ = new RSTFile(this->info_+this->path_+this->dir_,this->filename_);
	this->data_write_->precision(10);

	IOFiles corr_file(this->analyse_+this->path_+this->dir_+this->filename_+"-corr.dat",true);
	IOFiles long_range_corr_file(this->analyse_+this->path_+this->dir_+this->filename_+"-long-range-corr.dat",true);//should not be delcared when type!=2

	Vector<double> lrc_mean;
	unsigned int nruns;
	unsigned int tmax;

	(*this->read_)>>nruns>>tmax;
	(*this->data_write_)<<"% E dE 0|1"<<IOFiles::endl;
	/* the +1 is the averages over all runs */
	for(unsigned int i(0);i<nruns+1;i++){ 
		(*this->read_)>>this->E_>>this->corr_>>this->long_range_corr_;
		if(lrc_mean.size() == 0){ lrc_mean.set(this->long_range_corr_.size(),0);}
		(*this->data_write_)<<" "<<this->E_.get_x()<<" "<<this->E_.get_dx()<<" "<<(i<nruns?true:false)<<IOFiles::endl;
		for(unsigned int j(0);j<this->corr_.size();j++){
			corr_file<<j+0.5<<" "<<this->corr_[j]<<" "<<(i<nruns?true:false)<<IOFiles::endl;
		}
		for(unsigned int j(0);j<this->long_range_corr_.size();j++){
			long_range_corr_file<<j<<" "<<this->long_range_corr_[j]<<" "<<(i<nruns?true:false)<<IOFiles::endl;
			if(i<nruns){ lrc_mean(j) += this->long_range_corr_[j].get_x();}
		}
		if(this->long_range_corr_.size()>0){
			long_range_corr_file<<this->n_<<" "<<this->long_range_corr_[0]<<" "<<(i<nruns?true:false)<<IOFiles::endl;
		}
	}
	(*this->jd_write_)("E",this->E_);

	/*{*/
	Gnuplot gp(this->analyse_+this->path_+this->dir_,this->filename_+"-corr");
	gp+="set xlabel 'site' offset 0,0.5";
	gp+="set ylabel '$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(i+1)>$' offset 1";
	gp+="set title '$N="+tostring(this->N_)+"$ $m="+tostring(this->m_)+"$ $n="+tostring(this->n_)+"$ bc="+tostring(this->bc_)+"'";
	gp+="plot '"+this->filename_+"-corr.dat' u 1:($6==1?$2:1/0):3 w errorbars lt 1 lc 1 lw 2 t 'Independant measures',\\";
	gp+="     '"+this->filename_+"-corr.dat' u 1:($6==0?$2:1/0):3 w errorbars lt 1 lc 2 lw 2 t 'Mean'";
	gp.save_file();
	this->rst_file_->link_figure(this->analyse_+this->path_+this->dir_+this->filename_+"-corr.png","Correlation on links",this->analyse_+this->path_+this->dir_+this->filename_+"-corr.gp",1000);

	if(this->long_range_corr_.size()>0){
		unsigned int llr(this->long_range_corr_.size());
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

		IOFiles data_sf(this->analyse_+this->path_+this->dir_+this->filename_+"-structure-factor.dat",true);
		for(unsigned int k(0);k<llr;k++){
			data_sf<<dk*k<<" "<<Ck(k).real()<<" "<<Ck(k).imag()<<" "<<std::abs(Ck(k))<<IOFiles::endl;
		}

		Gnuplot gplr(this->analyse_+this->path_+this->dir_,this->filename_+"-long-range-corr");
		gplr.xrange(0,llr);
		gplr+="set xlabel '$\\|i-j\\|$' offset 0,0.5";
		gplr+="set ylabel '$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(j)>$' offset 1";
		gplr+="set title '$N="+tostring(this->N_)+"$ $m="+tostring(this->m_)+"$ $n="+tostring(this->n_)+"$ bc="+tostring(this->bc_)+"'";
		gplr+="set key right bottom";
		gplr+="a=1.0";
		gplr+="b=1.0";
		gplr+="eta=1.0";
		gplr+="m="+tostring(this->m_)+".0";
		gplr+="N="+tostring(this->N_)+".0";
		gplr+="plot '"+this->filename_+"-long-range-corr.dat' u 1:($6==1?$2:1/0):3 w errorbars lt 1 lc 1 lw 2 t 'Independant measures',\\";
		gplr+="     '"+this->filename_+"-long-range-corr.dat' u 1:($6==0?$2:1/0):3 w errorbars lt 1 lc 2 lw 2 t 'Mean'";
		gplr.save_file();
		gplr.create_image(true);
		this->rst_file_->link_figure(this->analyse_+this->path_+this->dir_+this->filename_+"-long-range-corr.png","Long range correlation",this->analyse_+this->path_+this->dir_+this->filename_+"-long-range-corr.gp",1000);

		Gnuplot gpsf(this->analyse_+this->path_+this->dir_,this->filename_+"-structure-factor");
		gpsf+="set title '$N="+tostring(this->N_)+"$ $m="+tostring(this->m_)+"$ $n="+tostring(this->n_)+"$ bc="+tostring(this->bc_)+"'";
		gpsf+="set key bottom";
		gpsf.xrange("0","2*pi");
		gpsf+="set xtics ('0' 0,'$\\pi/2$' pi/2,'$\\pi$' pi,'$3\\pi/2$' 3*pi/2,'$2\\pi$' 2 *pi)";
		gpsf+="plot '"+this->filename_+"-structure-factor.dat' u 1:2 t 'real',\\";
		gpsf+="     '"+this->filename_+"-structure-factor.dat' u 1:3 t 'imag',\\";
		gpsf+="     '"+this->filename_+"-structure-factor.dat' u 1:4 t 'norm'";
		gpsf.save_file();
		gpsf.create_image(true);
		this->rst_file_->link_figure(this->analyse_+this->path_+this->dir_+this->filename_+"-structure-factor.png","Structure factor",this->analyse_+this->path_+this->dir_+this->filename_+"-structure-factor.gp",1000);
	}
	/*}*/

	this->rst_file_->text(this->read_->get_header());
	this->rst_file_->save(false);
	delete this->rst_file_;
	this->rst_file_ = NULL;

	return this->filename_;
}

template<typename Type>
std::string ChainFermi<Type>::extract_level_6(){
	Data<double> min_E;
	min_E.set_x(0.0);

	unsigned int nof(0);
	(*this->read_)>>nof;
	for(unsigned int i(0);i<nof;i++){
		(*this->read_)>>this->E_;
		if(this->E_.get_x()<min_E.get_x()){ 
			min_E = this->E_;
		}
	}

	this->save();
	(*this->jd_write_)("energy per site",min_E);

	return this->filename_;
}

template<typename Type>
std::string ChainFermi<Type>::extract_level_3(){
	(*this->read_)>>this->E_;
	(*this->data_write_)<<this->n_<<" "<<this->E_<<IOFiles::endl;

	return this->filename_;
}
/*}*/
#endif
