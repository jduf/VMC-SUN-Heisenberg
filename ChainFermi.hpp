#ifndef DEF_CHAINFERMI
#define DEF_CHAINFERMI

/*{Creates a chain with uniform hopping parameter*/
/*!
 * To properly solve the degeneracy problem, this wavefunction selects the
 * eigenvector |E_F,-> or |E_F,->. This seems to have exactly the same effect
 * as choosing anti-periodic boundary conditions
 */
/*}*/
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
};

template<typename Type>
ChainFermi<Type>::ChainFermi(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc):
	System(ref,N,m,n,M,bc),
	Chain<Type>(1,"chain-fermi")
{
	if(this->status_==2){
		this->init_fermionic();

		this->system_info_.text("+ Spin chain with real and identical hopping");
		this->system_info_.text("  term between each sites :");
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
	this->compute_H();
	this->plot_band_structure();
	this->status_++;
}
/*}*/

/*{method needed for analysing*/
template<typename Type>
std::string ChainFermi<Type>::extract_level_7(){
	this->rst_file_ = new RSTFile(this->info_+this->path_+this->dir_,this->filename_);

	/*!extract jdbin*/
	/*{*/
	IOFiles corr_file(this->analyse_+this->path_+this->dir_+this->filename_+"-corr.dat",true);
	IOFiles lr_corr_file(this->analyse_+this->path_+this->dir_+this->filename_+"-long-range-corr.dat",true);

	Vector<double> lrc_mean(this->links_.row(),0);
	unsigned int nruns;
	unsigned int tmax;

	(*this->read_)>>nruns>>tmax;
	(*this->data_write_)<<"%E dx conv(0|1) #conv mean(0|1)"<<IOFiles::endl;
	corr_file<<"%(2i+1)/2 corr(i,i+1) dx conv(0|1) #conv mean(0|1)"<<IOFiles::endl;
	lr_corr_file<<"%j corr(i,j) dx conv(0|1) #conv mean(0|1)"<<IOFiles::endl;
	/*!the +1 is the average over all runs */
	for(unsigned int i(0);i<nruns+1;i++){ 
		(*this->read_)>>this->E_>>this->corr_>>this->lr_corr_;
		(*this->data_write_)<<this->E_<<" "<<(i<nruns)<<IOFiles::endl;
		for(unsigned int j(0);j<this->corr_.size();j++){
			corr_file<<j+0.5<<" "<<this->corr_[j]<<" "<<(i<nruns)<<IOFiles::endl;
		}
		for(unsigned int j(0);j<this->lr_corr_.size();j++){
			lr_corr_file<<j<<" "<<this->lr_corr_[j]<<" "<<(i<nruns)<<IOFiles::endl;
		}
		if(i<nruns){
			for(unsigned int j(0);j<this->lr_corr_.size();j++){
				lrc_mean(j) += this->lr_corr_[j].get_x()/nruns;
			}
		} else {
			for(unsigned int j(0);j<this->lr_corr_.size();j++){
				if(this->lr_corr_[j].get_conv()){ lrc_mean(j) = this->lr_corr_[j].get_x(); } 
			}
		}
	}
	/*}*/
	/*!nearest neighbourg correlations*/
	/*{*/
	Gnuplot gp(this->analyse_+this->path_+this->dir_,this->filename_+"-corr");
	gp+="set xlabel 'site' offset 0,0.5";
	gp+="set y2label '$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(i+1)>$'";
	gp+="set title '$N="+tostring(this->N_)+"$ $m="+tostring(this->m_)+"$ $n="+tostring(this->n_)+"$ bc="+tostring(this->bc_)+"'";
	gp+="plot '"+this->filename_+"-corr.dat' u 1:(($6==1 && $5==0)?$2:1/0):3 w errorbars lt 1 lc 5 t 'Not converged',\\";
	gp+="     '"+this->filename_+"-corr.dat' u 1:(($6==1 && $5==1)?$2:1/0):3 w errorbars lt 1 lc 6 t 'Converged',\\";
	gp+="     '"+this->filename_+"-corr.dat' u 1:($6==0?$2:1/0):3 w errorbars lt 1 lc 7 t 'Mean'";
	gp.save_file();
	//gp.create_image(true);
	this->rst_file_->link_figure(this->analyse_+this->path_+this->dir_+this->filename_+"-corr.png","Correlation on links",this->analyse_+this->path_+this->dir_+this->filename_+"-corr.gp",1000);
	/*}*/
	/*!long range correlations*/
	/*{*/
	unsigned int xi(1);
	unsigned int xf(this->n_);
	Vector<double> exponents(4);
	this->compute_critical_exponents(xi,xf,exponents,lrc_mean);

	Gnuplot gplr(this->analyse_+this->path_+this->dir_,this->filename_+"-long-range-corr");
	gplr.xrange(this->N_/this->m_,this->n_-this->N_/this->m_);
	gplr+="set xlabel '$\\|i-j\\|$' offset 0,0.5";
	gplr+="set y2label '$<S_{\\alpha}^{\\alpha}(i)S_{\\alpha}^{\\alpha}(j)>-\\dfrac{m^2}{N}$' offset 1";
	gplr+="set title '$N="+tostring(this->N_)+"$ $m="+tostring(this->m_)+"$ $n="+tostring(this->n_)+"$ bc="+tostring(this->bc_)+"'";
	gplr+="set key center bottom";
	gplr+="set sample 1000";
	gplr+="m="+tostring(this->m_)+".0";
	gplr+="N="+tostring(this->N_)+".0";
	gplr+="n="+tostring(this->n_)+".0";
	gplr+="p0 = 1.0";
	gplr+="p1 = 2.0-2.0/N";
	gplr+="p2 = -1.0";
	gplr+="p3 = 2.0";
	gplr+="f(x) = p0*cos(2.0*pi*x*m/N)*(x**(-p1)+(n-x)**(-p1))+p2*(x**(-p3)+(n-x)**(-p3))";
	gplr+="set fit quiet";
	gplr+="fit [" + tostring(xi) + ":" + tostring(xf) + "] f(x) '"+this->filename_+"-long-range-corr.dat' u 1:($6==0?$2:1/0) noerrors via p0,p1,p2,p3"; 
	gplr+="plot '"+this->filename_+"-long-range-corr.dat' u 1:(($6==1 && $5==0)?$2:1/0):3 w errorbars lt 1 lc 5 t 'Not converged',\\";
	gplr+="     '"+this->filename_+"-long-range-corr.dat' u 1:(($6==1 && $5==1)?$2:1/0):3 w errorbars lt 1 lc 6 t 'Converged',\\";
	gplr+="     '"+this->filename_+"-long-range-corr.dat' u 1:($6==0?$2:1/0):3 w errorbars lt 1 lc 7 t 'Mean',\\";
	gplr+="     f(x) lc 7 lw 0.5 t sprintf('$\\eta=%f$, $\\mu=%f$',p1,p3)";
	gplr.save_file();
//	gplr.create_image(true);
	this->rst_file_->link_figure(this->analyse_+this->path_+this->dir_+this->filename_+"-long-range-corr.png","Long range correlation",this->analyse_+this->path_+this->dir_+this->filename_+"-long-range-corr.gp",1000);
	/*}*/
	/*!structure factor*/
	/*{*/
	unsigned int llr(this->lr_corr_.size());
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

	IOFiles data_sf(this->analyse_+this->path_+this->dir_+this->filename_+"-structure-factor.dat",true);
	for(unsigned int k(0);k<llr;k++){
		data_sf<<dk*k<<" "<<Ck(k).real()<<" "<<Ck(k).imag()<<IOFiles::endl;
	}

	Gnuplot gpsf(this->analyse_+this->path_+this->dir_,this->filename_+"-structure-factor");
	gpsf+="set title '$N="+tostring(this->N_)+"$ $m="+tostring(this->m_)+"$ $n="+tostring(this->n_)+"$ bc="+tostring(this->bc_)+"'";
	gpsf+="set key bottom";
	gpsf.xrange("0","2*pi");
	switch(this->N_/this->m_){
		case 3: { gpsf+="set xtics ('0' 0,'$2\\pi/3$' 2.0*pi/3.0, '$4\\pi/3$' 4.0*pi/3.0,'$2\\pi$' 2.0*pi)"; } break;
		case 5: { gpsf+="set xtics ('0' 0,'$2\\pi/5$' 2.0*pi/5.0, '$4\\pi/5$' 4.0*pi/5.0, '$6\\pi/5$' 6.0*pi/5.0, '$8\\pi/5$' 8.0*pi/5.0, '$2\\pi$' 2.0*pi)"; } break;
		default:{ gpsf+="set xtics ('0' 0,'$\\pi/2$' pi/2.0,'$\\pi$' pi,'$3\\pi/2$' 3.0*pi/2.0,'$2\\pi$' 2.0*pi)"; } break;
	}
	gpsf+="set xlabel '$k$' offset 0,0.5";
	gpsf+="set y2label '$<S(k)>$'";
	gpsf+="plot '"+this->filename_+"-structure-factor.dat' u 1:2 lt 1 lc 6 t 'real',\\";
	gpsf+="     '"+this->filename_+"-structure-factor.dat' u 1:3 lt 1 lc 7 t 'imag'";
	gpsf.save_file();
	//gpsf.create_image(true);
	this->rst_file_->link_figure(this->analyse_+this->path_+this->dir_+this->filename_+"-structure-factor.png","Structure factor",this->analyse_+this->path_+this->dir_+this->filename_+"-structure-factor.gp",1000);
	/*}*/
	/*!save some additionnal values */
	/*{*/
	this->jd_write_->write("energy per site",this->E_);
	this->jd_write_->write("polymerization strength",0.0);
	this->jd_write_->write("critical exponents",exponents);
	/*}*/

	this->rst_file_->text(this->read_->get_header());
	this->rst_file_->save(false);
	delete this->rst_file_;
	this->rst_file_ = NULL;

	return this->filename_;
}

template<typename Type>
std::string ChainFermi<Type>::extract_level_6(){
	double polymerization_strength;
	Vector<double> exponents;
	unsigned int nof(0);
	(*this->read_)>>nof>>this->E_>>polymerization_strength>>exponents;

	this->jd_write_->add_to_header("\n");
	this->save();
	this->jd_write_->write("energy per site",this->E_);
	this->jd_write_->write("polymerization strength",polymerization_strength);
	this->jd_write_->write("critical exponents",exponents);

	return this->filename_;
}
/*}*/
#endif
