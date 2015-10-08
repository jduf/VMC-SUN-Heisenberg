#ifndef DEF_CHAINFERMI
#define DEF_CHAINFERMI

#include "Chain.hpp"

/*{Description*/
/*!Creates a chain with uniform hopping parameter
 * To properly solve the degeneracy problem, this wavefunction selects the
 * eigenvector |E_F,-> or |E_F,->. This seems to have exactly the same effect
 * as choosing anti-periodic boundary conditions
 */
/*}*/
template<typename Type>
class ChainFermi: public Chain<Type>{
	public:
		ChainFermi(System const& s);
		~ChainFermi() = default;

		void create();
		void check();

	private:
		void compute_H();
		std::string extract_level_8();
		std::string extract_level_7();
};

template<typename Type>
ChainFermi<Type>::ChainFermi(System const& s):
	System(s),
	Chain<Type>(1,"chain-fermi")
{
	if(this->status_==2){
		this->init_fermionic();

		this->system_info_.item("+ Spin chain with real and identical hopping term between each sites");
	}
}

/*{method needed for running*/
template<typename Type>
void ChainFermi<Type>::compute_H(){
	this->H_.set(this->n_,this->n_,0);
	Matrix<int> nb;
	for(unsigned int i(0);i<this->n_;i++){
		nb = this->get_neighbourg(i);
		this->H_(i,nb(0,0)) = nb(0,1);
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
std::string ChainFermi<Type>::extract_level_8(){
	this->rst_file_ = new RSTFile(this->info_+this->path_+this->dir_,this->filename_);
	std::string basename("../../../../../../../../"+this->analyse_+this->path_+this->dir_+this->filename_);
	std::string title("$N="+my::tostring(this->N_)+"$ $m="+my::tostring(this->m_)+"$ $n="+my::tostring(this->n_)+"$ bc="+my::tostring(this->bc_));

	/*!extract jdbin*/
	/*{*/
	IOFiles corr_file(this->analyse_+this->path_+this->dir_+this->filename_+"-corr.dat",true);
	IOFiles lr_corr_file(this->analyse_+this->path_+this->dir_+this->filename_+"-long-range-corr.dat",true);

	Vector<double> lr_corr_v(this->links_.row());

	corr_file<<"%(2i+1)/2 corr(i,i+1) dx conv(0|1) #conv mean(0|1)"<<IOFiles::endl;
	lr_corr_file<<"%j corr(i,j) dx conv(0|1) #conv mean(0|1)"<<IOFiles::endl;

	(*this->read_)>>this->E_>>this->corr_>>this->lr_corr_;
	(*this->data_write_)<<this->E_<<IOFiles::endl;
	for(unsigned int i(0);i<this->corr_.size();i++){
		corr_file<<i+0.5<<" "<<this->corr_[i]<<IOFiles::endl;
	}
	for(unsigned int i(0);i<this->lr_corr_.size();i++){
		lr_corr_file<<i<<" "<<this->lr_corr_[i]<<IOFiles::endl;
		lr_corr_v(i) = this->lr_corr_[i].get_x();
	}
	/*}*/
	/*!nearest neighbourg correlations*/
	/*{*/
	Gnuplot gp(this->analyse_+this->path_+this->dir_,this->filename_+"-corr");
	gp.label("x","site","offset 0,0.5");
	gp.label("y2","$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(i+1)>$");
	gp.title(title);
	gp+="plot '"+this->filename_+"-corr.dat' u 1:2:3 w errorbars lt 1 lc 7 notitle";
	gp.save_file();
	//gp.create_image(true);
	this->rst_file_->figure(basename+"-corr.png","Correlation on links",RST::target(basename+"-corr.gp")+RST::width("1000"));
	/*}*/
	/*!long range correlations*/
	/*{*/
	unsigned int xi;
	unsigned int xf;
	Vector<double> exponents;
	bool fit(this->compute_critical_exponents(lr_corr_v,xi,xf,exponents));

	Gnuplot gplr(this->analyse_+this->path_+this->dir_,this->filename_+"-long-range-corr");
	gplr.range("x",this->N_/this->m_,this->n_-this->N_/this->m_);
	gplr.label("x","$\\|i-j\\|$","offset 0,0.5");
	gplr.label("y2","$<S_{\\alpha}^{\\alpha}(i)S_{\\alpha}^{\\alpha}(j)>-\\dfrac{m^2}{N}$","offset 1");
	gp.title(title);
	gplr+="set key center bottom";
	gplr+="set sample 1000";
	gplr+="m="+my::tostring(this->m_)+".0";
	gplr+="N="+my::tostring(this->N_)+".0";
	gplr+="n="+my::tostring(this->n_)+".0";
	gplr+="p0 = 1.0";
	gplr+="p1 = 2.0-2.0/N";
	gplr+="p2 = -1.0";
	gplr+="p3 = 2.0";
	gplr+="f(x) = p0*cos(2.0*pi*x*m/N)*(x**(-p1)+(n-x)**(-p1))+p2*(x**(-p3)+(n-x)**(-p3))";
	gplr+="set fit quiet";
	gplr+="fit [" + my::tostring(xi) + ":" + my::tostring(xf) + "] f(x) '"+this->filename_+"-long-range-corr.dat' u 1:($6==0?$2:1/0) noerrors via p0,p1,p2,p3"; 
	gplr+="plot '"+this->filename_+"-long-range-corr.dat' u 1:2:3 w errorbars lt 1 lc 7 notitle,\\";
	gplr+="     f(x) lc 7 " + std::string(fit?"lw 0.5":"dt 2") + " t sprintf('$\\eta=%f$, $\\mu=%f$',p1,p3)";
	gplr.save_file();
	//gplr.create_image(true);
	this->rst_file_->figure(basename+"-long-range-corr.png","Long range correlation",RST::target(basename+"-long-range-corr.gp")+RST::width("1000"));
	/*}*/
	/*!structure factor*/
	/*{*/
	unsigned int llr(this->lr_corr_.size());
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

	IOFiles data_sf(this->analyse_+this->path_+this->dir_+this->filename_+"-structure-factor.dat",true);
	for(unsigned int k(0);k<llr;k++){
		data_sf<<dk*k<<" "<<Ck(k).real()<<" "<<Ck(k).imag()<<IOFiles::endl;
	}

	Gnuplot gpsf(this->analyse_+this->path_+this->dir_,this->filename_+"-structure-factor");
	gp.title(title);
	gpsf+="set key bottom";
	gpsf.range("x","0","2*pi");
	switch(this->N_/this->m_){
		case 3: { gpsf+="set xtics ('0' 0,'$2\\pi/3$' 2.0*pi/3.0, '$4\\pi/3$' 4.0*pi/3.0,'$2\\pi$' 2.0*pi)"; } break;
		case 5: { gpsf+="set xtics ('0' 0,'$2\\pi/5$' 2.0*pi/5.0, '$4\\pi/5$' 4.0*pi/5.0, '$6\\pi/5$' 6.0*pi/5.0, '$8\\pi/5$' 8.0*pi/5.0, '$2\\pi$' 2.0*pi)"; } break;
		default:{ gpsf+="set xtics ('0' 0,'$\\pi/2$' pi/2.0,'$\\pi$' pi,'$3\\pi/2$' 3.0*pi/2.0,'$2\\pi$' 2.0*pi)"; } break;
	}
	gpsf.label("x","$k$","offset 0,0.5");
	gpsf.label("y2","$<S(k)>$");
	gpsf+="plot '"+this->filename_+"-structure-factor.dat' u 1:2 lt 1 lc 6 t 'real',\\";
	gpsf+="     '"+this->filename_+"-structure-factor.dat' u 1:3 lt 1 lc 7 t 'imag'";
	gpsf.save_file();
	//gpsf.create_image(true);
	this->rst_file_->figure(basename+"-structure-factor.png","Structure factor",RST::target(basename+"-structure-factor.gp")+RST::width("1000"));
	/*}*/
	/*!save*/
	/*{*/
	this->jd_write_->add_header()->title("System's parameters",'-');
	this->save_param(*this->jd_write_);
	this->save_input(*this->jd_write_);
	this->save_output(*this->jd_write_);
	this->jd_write_->write("polymerization strength",0.0);
	this->jd_write_->write("critical exponents",exponents);
	/*}*/

	this->rst_file_->text(this->read_->get_header());
	this->rst_file_->save(false,true);
	delete this->rst_file_;
	this->rst_file_ = NULL;

	return this->filename_;
}

template<typename Type>
std::string ChainFermi<Type>::extract_level_7(){

	this->jd_write_->add_header()->title("System's parameters",'-');
	this->save_param(*this->jd_write_);
	this->save_input(*this->jd_write_);
	this->save_output(*this->jd_write_);
	this->jd_write_->write("energy per site",this->E_);
	double polymerization_strength;
	Vector<double> exponents;
	(*this->read_)>>polymerization_strength>>exponents;
	this->jd_write_->write("polymerization strength",polymerization_strength);
	this->jd_write_->write("critical exponents",exponents);

	return this->filename_;
}
/*}*/
#endif
