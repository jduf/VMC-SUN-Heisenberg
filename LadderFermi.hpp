#ifndef DEF_LADDERFERMI
#define DEF_LADDERFERMI

#include "Ladder.hpp"

/*{Description*/
/*!Creates a ladder with uniform hopping parameter 
 * 
 *  => Fermi Ladder<=
 * 
 * */
/*}*/
template<typename Type>
class LadderFermi: public Ladder<Type>{
	public:
		LadderFermi(System const& s);
		~LadderFermi(){}

		void create();
		void check();

	private:
		void compute_H();
		void lattice(std::string const& path, std::string const& filename);

		std::string extract_level_7();
		std::string extract_level_6();
};

template<typename Type>
LadderFermi<Type>::LadderFermi(System const& s):
	System(s),
	Ladder<Type>(2,"ladder-fermi")
{
	if(this->status_==2){
		this->init_fermionic();

		this->system_info_.item("Spin ladder with real and identical hopping term between each sites :");
		this->system_info_.item("   => Fermi ladder <=   ");
	}
}

/*{method needed for running*/
template<typename Type>
void LadderFermi<Type>::compute_H(){
	this->H_.set(this->n_,this->n_,0);
	Matrix<int> nb;
	for(unsigned int i(0);i<this->n_;i++){
		nb = this->get_neighbourg(i);
		if( ( i % 2 ) == 0){				//even sites => top
			this->H_(i,nb(0,0)) = nb(0,1);
			this->H_(i,nb(1,0)) = nb(1,1);	
		} else {							//odd sites => bottom
			this->H_(i,nb(0,0)) = nb(0,1);
		}
	}
	this->H_ += this->H_.transpose();
}
/*}*/

/*{method needed for checking*/
template<typename Type>
void LadderFermi<Type>::check(){
	this->compute_H();

	Matrix<int> nb;
	for(unsigned int i(0);i<this->n_;i++){
		nb = this->get_neighbourg(i);
		std::cout<<"i="<<i;
		std::cout<<" nb="<<nb(0,0);
		std::cout<<","<<nb(1,0);
		std::cout<<","<<nb(2,0)<<std::endl;;
	} 
	for(unsigned int i(0);i<this->links_.row();i++){
		std::cout<<"i="<<i;
		std::cout<<" l="<<this->links_(i,0);
		std::cout<<","<<this->links_(i,1);
		std::cout<<" J="<<this->J_(i)<<std::endl;;
	} 
	std::cout<<"Hamiltonien"<<std::endl;
	std::cout<<this->H_<<std::endl;

	this->plot_band_structure();
	lattice("./","lattice");
}

template<typename Type>
void LadderFermi<Type>::lattice(std::string const& path, std::string const& filename){
	compute_H();
	Matrix<int> nb;
	std::string color("black");
	std::string linestyle("solid");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(path,filename);
	ps.begin(-9,-10,16,10,this->filename_);
	for(unsigned int i(0);i<this->n_;i++) {
		xy0(0) = i/2;
		xy0(1) = i%2;
		ps.put(xy0(0)-0.20,xy0(1)+0.15,my::tostring(i));
		nb = this->get_neighbourg(i);

		if(nb(0,1)<0){ color = "red"; } 
		else { color = "black"; }
		xy1(0) = nb(0,0)/2;
		xy1(1) = nb(0,0)%2;
		if(xy1(0)<xy0(0)){ 
			xy1(0) = xy0(0)+1;
			linestyle="dashed";
		} else{ linestyle="solid"; }
		/*x-link*/  ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linecolor="+color+",linestyle="+linestyle);
		if(i%2){
			color = "black";
			linestyle="solid"; 
			xy1(0) = nb(1,0)/2;
			xy1(1) = nb(1,0)%2;
			/*y-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linecolor="+color+",linestyle="+linestyle);
		}
	}

	Matrix<double> polygon(4,2);
	polygon(0,0)=-0.1;
	polygon(0,1)=-0.1;
	polygon(1,0)=this->n_/2-0.1;
	polygon(1,1)=-0.1;
	polygon(2,0)=this->n_/2-0.1;
	polygon(2,1)=1.1;
	polygon(3,0)=-0.1;
	polygon(3,1)=1.1;
	ps.polygon(polygon,"linecolor=green");

	polygon(0,0)=-0.1;
	polygon(0,1)=-0.1;
	polygon(1,0)=0.9;
	polygon(1,1)=-0.1;
	polygon(2,0)=0.9;
	polygon(2,1)=1.1;
	polygon(3,0)=-0.1;
	polygon(3,1)=1.1;
	ps.polygon(polygon,"linecolor=blue");

	ps.end(true,true,true);
}
/*}*/

/*{method needed for analysing*/
template<typename Type>
std::string LadderFermi<Type>::extract_level_7(){
	this->rst_file_ = new RSTFile(this->info_+this->path_+this->dir_,this->filename_);
	std::string title("$N="+my::tostring(this->N_)+"$ $m="+my::tostring(this->m_)+"$ $n="+my::tostring(this->n_)+"$ bc="+my::tostring(this->bc_));

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
	gp.label("x","site","offset 0,0.5");
	gp.label("y2","$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(i+1)>$");
	gp.title(title);
	gp+="plot '"+this->filename_+"-corr.dat' u 1:(($6==1 && $5==0)?$2:1/0):3 w errorbars lt 1 lc 5 t 'Not converged',\\";
	gp+="     '"+this->filename_+"-corr.dat' u 1:(($6==1 && $5==1)?$2:1/0):3 w errorbars lt 1 lc 6 t 'Converged',\\";
	gp+="     '"+this->filename_+"-corr.dat' u 1:($6==0?$2:1/0):3 w errorbars lt 1 lc 7 t 'Mean'";
	gp.save_file();
	//gp.create_image(true);
	this->rst_file_->figure(this->analyse_+this->path_+this->dir_+this->filename_+"-corr.png","Correlation on links",RST::target(this->analyse_+this->path_+this->dir_+this->filename_+"-corr.gp")+RST::width("1000"));
	/*}*/
	/*!long range correlations*/
	/*{*/
	Vector<double> exponents;

	Gnuplot gplr(this->analyse_+this->path_+this->dir_,this->filename_+"-long-range-corr");
	gplr.range("x",this->N_/this->m_,this->n_-this->N_/this->m_);
	gplr.label("x","$\\|i-j\\|$","offset 0,0.5");
	gplr.label("y2","$<S_{\\alpha}^{\\alpha}(i)S_{\\alpha}^{\\alpha}(j)>-\\dfrac{m^2}{N}$","offset 1");
	gp.title(title);
	gplr+="set key center bottom";
	gplr+="m="+my::tostring(this->m_)+".0";
	gplr+="N="+my::tostring(this->N_)+".0";
	gplr+="n="+my::tostring(this->n_)+".0";
	gplr+="plot '"+this->filename_+"-long-range-corr.dat' u 1:(($6==1 && $5==0)?$2:1/0):3 w errorbars lt 1 lc 5 t 'Not converged',\\";
	gplr+="     '"+this->filename_+"-long-range-corr.dat' u 1:(($6==1 && $5==1)?$2:1/0):3 w errorbars lt 1 lc 6 t 'Converged',\\";
	gplr+="     '"+this->filename_+"-long-range-corr.dat' u 1:($6==0?$2:1/0):3 w errorbars lt 1 lc 7 t 'Mean'";
	gplr.save_file();
	//gplr.create_image(true);
	this->rst_file_->figure(this->analyse_+this->path_+this->dir_+this->filename_+"-long-range-corr.png","Long range correlation",RST::target(this->analyse_+this->path_+this->dir_+this->filename_+"-long-range-corr.gp")+RST::width("1000"));
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
	gp.title(title);
	gpsf+="set key bottom";
	gpsf.range("x","0","2*pi");
	gpsf+="set xtics ('0' 0,'$\\pi/2$' pi/2.0,'$\\pi$' pi,'$3\\pi/2$' 3.0*pi/2.0,'$2\\pi$' 2.0*pi)"; 
	gpsf.label("x","$k$","offset 0,0.5");
	gpsf.label("y2","$<S(k)>$");
	gpsf+="plot '"+this->filename_+"-structure-factor.dat' u 1:2 lt 1 lc 6 t 'real',\\";
	gpsf+="     '"+this->filename_+"-structure-factor.dat' u 1:3 lt 1 lc 7 t 'imag'";
	gpsf.save_file();
	//gpsf.create_image(true);
	this->rst_file_->figure(this->analyse_+this->path_+this->dir_+this->filename_+"-structure-factor.png","Structure factor",RST::target(this->analyse_+this->path_+this->dir_+this->filename_+"-structure-factor.gp")+RST::width("1000"));
	/*}*/
	/*!save some additionnal values */
	/*{*/
	this->jd_write_->write("energy per site",this->E_);
	/*}*/

	this->rst_file_->text(this->read_->get_header());
	this->rst_file_->save(false,true);
	delete this->rst_file_;
	this->rst_file_ = NULL;

	return this->filename_;
}

template<typename Type>
std::string LadderFermi<Type>::extract_level_6(){
	unsigned int nof(0);
	(*this->read_)>>nof>>this->E_;

	this->jd_write_->add_header()->nl();
	this->save_input(*this->jd_write_);
	this->jd_write_->write("energy per site",this->E_);

	return this->filename_;
}
/*}*/
#endif
