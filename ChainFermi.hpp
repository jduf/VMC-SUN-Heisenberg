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
		std::string extract_level_4();
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
	create();
	//BandStructure<Type> bs(this->compute_T(),this->Lx_,this->spuc_,this->bc_);
}
/*}*/

/*{method needed for analysing*/
template<typename Type>
std::string ChainFermi<Type>::extract_level_7(){
	this->rst_file_ = new RSTFile(this->info_+this->path_+this->dir_,this->filename_);

	unsigned int nruns;
	unsigned int tmax;

	(*this->read_)>>nruns>>tmax;
	//std::cout<<delta_<<" "<<this->N_<<" "<<this->m_<<" "<<this->n_<<" "<<this->M_<<" "<<tmax<<" "<<nruns<<std::endl;
	IOFiles corr_file(this->analyse_+this->path_+this->dir_+this->filename_+"-corr.dat",true);
	IOFiles long_range_corr_file(this->analyse_+this->path_+this->dir_+this->filename_+"-long-range-corr.dat",true);//should not be delcared when type!=2
	this->data_write_->precision(10);
	(*this->data_write_)<<"% E dE 0|1"<<IOFiles::endl;
	/* the +1 is the averages over all runs */
	Vector<double> poly_e(this->N_/this->m_,0);
	for(unsigned int i(0);i<nruns+1;i++){ 
		(*this->read_)>>this->E_>>this->corr_>>this->long_range_corr_;
		if(i<nruns){
			unsigned int k(0);
			while(k<this->corr_.size()){
				for(unsigned int j(0);j<this->N_/this->m_;j++){
					poly_e(j) += this->corr_[k].get_x();
					k++;
				}
			}
		}
		(*this->data_write_)<<" "<<this->E_.get_x()<<" "<<this->E_.get_dx()<<" "<<(i<nruns?true:false)<<IOFiles::endl;
		for(unsigned int j(0);j<this->corr_.size();j++){
			corr_file<<j+0.5<<" "<<this->corr_[j]<<" "<<(i<nruns?true:false)<<IOFiles::endl;
		}
		for(unsigned int j(0);j<this->long_range_corr_.size();j++){
			long_range_corr_file<<j+1<<" "<<this->long_range_corr_[j]<<" "<<(i<nruns?true:false)<<IOFiles::endl;
		}
	}
	poly_e /= nruns*this->n_*this->m_/this->N_;
	poly_e.sort(std::less<double>());

	(*this->jd_write_)("E",this->E_);
	(*this->jd_write_)("polymerization strength",poly_e(this->N_/this->m_-1)-poly_e(this->N_/this->m_-2));

	/*{*/
	Gnuplot gp(this->analyse_+this->path_+this->dir_,this->filename_+"-corr");
	gp+="set xlabel 'site' offset 0,0.5";
	gp+="set ylabel '$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(i+1)>$' offset 1";
	gp+="set title '$N="+tostring(this->N_)+"$ $m="+tostring(this->m_)+"$ $n="+tostring(this->n_)+"$ bc="+tostring(this->bc_)+"'";
	gp+="plot '"+this->filename_+"-corr.dat' u 1:($6==1?$2:1/0):3 w errorbars lt 1 lc 1 lw 2 t 'Independant measures',\\";
	gp+="     '"+this->filename_+"-corr.dat' u 1:($6==0?$2:1/0):3 w errorbars lt 1 lc 2 lw 2 t 'Mean',\\";
	gp+="     "+tostring(poly_e(this->N_/this->m_-1)) + " w l lc 3 t 'd-merization="+tostring(poly_e(this->N_/this->m_-1)-poly_e(this->N_/this->m_-2))+"',\\";
	gp+="     "+tostring(poly_e(this->N_/this->m_-2)) + " w l lc 3 notitle";
	gp.save_file();
	this->rst_file_->link_figure(this->analyse_+this->path_+this->dir_+this->filename_+"-corr.png","Correlation on links",this->analyse_+this->path_+this->dir_+this->filename_+"-corr.gp",1000);

	unsigned int length(this->long_range_corr_.size());
	if(length>0){
		Gnuplot gp(this->analyse_+this->path_+this->dir_,this->filename_+"-long-range-corr");
		gp+="stats '"+this->filename_+"-long-range-corr.dat' nooutput";
		gp.xrange(0,length+1);
		gp.yrange("1.1*STATS_min_y","1.1*STATS_max_y");
		gp+="set xlabel '$\\|i-j\\|$' offset 0,0.5";
		gp+="set ylabel '$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(j)>$' offset 1";
		gp+="set title '$N="+tostring(this->N_)+"$ $m="+tostring(this->m_)+"$ $n="+tostring(this->n_)+"$ bc="+tostring(this->bc_)+"'";
		gp+="set key right bottom";
		gp+="a=1.0";
		gp+="b=1.0";
		gp+="eta=1.0";
		gp+="m="+tostring(this->m_)+".0";
		gp+="N="+tostring(this->N_)+".0";
		gp+="f(x) = a/(x*x) + b*cos(2.0*pi*x*m/N)/(x**eta)";
		gp+="set fit quiet";
		switch(this->N_/this->m_){
			case 2:{ gp+="fit [3:"+tostring(length)+"] f(x) '"+this->filename_+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; } break;
			case 3:{
					   switch((length + 1) % 3){
						   case 0:{ gp+="fit [2:"+tostring(length)+"] f(x) '"+this->filename_+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
						   case 1:{ gp+="fit [5:"+tostring(length)+"] f(x) '"+this->filename_+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
						   case 2:{ gp+="fit [3:"+tostring(length)+"] f(x) '"+this->filename_+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
					   }break;
				   }break;
			default :{ gp+="fit ["+tostring(this->N_-1)+":] f(x) '"+this->filename_+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
		}
		gp+="plot for [IDX=0:"+tostring(nruns-1)+"] '"+this->filename_+"-long-range-corr.dat' i IDX u 1:($4==1?$2:1/0):3 w errorbars lt 1 lc 3 ps 0 notitle,\\";
		gp+="     for [IDX=0:"+tostring(nruns-1)+"] '"+this->filename_+"-long-range-corr.dat' i IDX u 1:($4==0?$2:1/0):3 w errorbars lt 1 lc 1 ps 0 notitle,\\";
		gp+="                   '"+this->filename_+"-long-range-corr.dat' i "+tostring(nruns)+" u 1:2:3 w errorbars lt 1 lc 2 lw 2 notitle,\\";
		gp+="                   f(x) notitle";
		gp.save_file();
		this->rst_file_->link_figure(this->analyse_+this->path_+this->dir_+this->filename_+"-long-range-corr.png","Long range correlation",this->analyse_+this->path_+this->dir_+this->filename_+"-long-range-corr.gp",1000);
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
	double min_polymerization_strength(0.0);
	double polymerization_strength;
	Data<double> min_E;
	min_E.set_x(0.0);

	unsigned int nof(0);
	(*this->read_)>>nof;
	for(unsigned int i(0);i<nof;i++){
		(*this->read_)>>this->E_>>polymerization_strength;
		if(this->E_.get_x()<min_E.get_x()){ 
			min_E = this->E_;
			min_polymerization_strength = polymerization_strength;
		}
	}

	this->save();
	(*this->jd_write_)("energy per site",min_E);
	(*this->jd_write_)("polymerization strength",min_polymerization_strength);

	return this->filename_;
}

template<typename Type>
std::string ChainFermi<Type>::extract_level_4(){
	double polymerization_strength;
	(*this->read_)>>this->E_>>polymerization_strength;

	this->save();
	(*this->jd_write_)("energy per site",this->E_);
	(*this->jd_write_)("polymerization strength",polymerization_strength);

	return this->filename_;
}

template<typename Type>
std::string ChainFermi<Type>::extract_level_3(){
	double polymerization_strength;
	(*this->read_)>>this->E_>>polymerization_strength;
	(*this->data_write_)<<this->n_<<" "<<polymerization_strength<<IOFiles::endl;

	return this->filename_;
}
/*}*/
#endif
