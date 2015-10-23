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

		void energy_bound(std::string const& path, std::string const& title);
};

template<typename Type>
ChainFermi<Type>::ChainFermi(System const& s):
	System(s),
	Chain<Type>(1,"chain-fermi")
{
	if(this->status_==2){
		this->init_fermionic();

		this->system_info_.item("Spin chain with real and identical hopping term between each sites");
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

template<typename Type>
void ChainFermi<Type>::energy_bound(std::string const& path, std::string const& filename){
	IOFiles corr_file(path+filename+"-corr.dat",true);
	corr_file<<"%(2i+1)/2 corr(i,i+1) dx conv(0|1) #conv mean(0|1)"<<IOFiles::endl;

	for(unsigned int i(0);i<this->obs_[0].nval();i++){
		corr_file<<i+0.5<<" "<<this->obs_[0][i]<<IOFiles::endl;
	}

	Gnuplot gp(path,filename+"-corr");
	gp+="set key center";
	gp.label("x","site","offset 0,0.5");
	gp.label("y2","$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(i+1)>$");
	gp+="plot '"+filename+"-corr.dat' u 1:2:3 w errorbars lt 1 lc 7 notitle";
	gp.save_file();
	gp.create_image(true,true);

	if(this->jd_write_){ this->jd_write_->write("polymerization strength",0.0); }
}
/*}*/

/*{method needed for analysing*/
template<typename Type>
std::string ChainFermi<Type>::extract_level_8(){
	this->rst_file_ = new RSTFile(this->info_+this->path_+this->dir_,this->filename_);
	std::string basename("../../../../../../../../"+this->analyse_+this->path_+this->dir_+this->filename_);
	std::string title("$N="+my::tostring(this->N_)+"$ $m="+my::tostring(this->m_)+"$ $n="+my::tostring(this->n_)+"$ bc="+my::tostring(this->bc_));

	(*this->data_write_)<<this->E_<<IOFiles::endl;
	this->jd_write_->add_header()->title("System's parameters",'-');
	this->save_param(*this->jd_write_);
	this->save_input(*this->jd_write_);
	this->save_output(*this->jd_write_);

	energy_bound(this->analyse_+this->path_+this->dir_,this->filename_);
	this->rst_file_->figure(basename+"-corr.png","Correlation on links",RST::target(basename+"-corr.gp")+RST::width("1000"));

	this->long_range_correlation_and_structure_factor(this->analyse_+this->path_+this->dir_,this->filename_);
	this->rst_file_->figure(basename+"-long-range-corr.png","Long range correlation",RST::target(basename+"-long-range-corr.gp")+RST::width("1000"));
	this->rst_file_->figure(basename+"-structure-factor.png","Structure factor",RST::target(basename+"-structure-factor.gp")+RST::width("1000"));

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
