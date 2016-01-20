#ifndef DEF_CHAIN
#define DEF_CHAIN

#include "System1D.hpp"

/*{Description*/
/*!1D chain with all hopping term having the same amplitude. The band structure
 * looks like this :
 *
 *     n_ even, bc_ = 1 |    bc_ = -1
 *                      |
 *           +          |      + +
 *         +   +        |    +     +
 *       +       +      |  +         +
 *                 + (1)|(3)
 *     -----------------|---------------
 *     n_ odd, bc_ = 1  |    bc_ = -1
 *                      |
 *           +          |      + +
 *         +   +        |    +     +
 *       +       +   (2)|(4)         +
 *     -----------------|---------------
 *
 * For the polymerized case, with N/m=spuc and spuc integer, unit cell contains
 * spuc sites and the brioullin zone is accordingly reduced. spuc!=1 when
 * di/tri/...-merization is created by ChainPolymerized with different hopping
 * term every spuc sites (delta!=0). It those cases, the selection of
 * eigenvector is unequivocal and there is no need to worry further. But when
 * delta==0, ChainFermi and ChainPolymerized are equivalent and may suffers
 * from the same problem, the selection of eigenvector may be equivocal because
 * at the Fermi level, the energies are degenerate :
 *
 * + N/m is odd, always degenerate but no discontinuity when delta->0 if the
 * |0> is selected
 * + N/m is even, everything works properly when nm/N is odd but if it is even,
 * then E(|0>) > E(delta->0) > E(|E_F,+>)
 */
/*}*/
template<typename Type>
class Chain: public System1D<Type>{
	public:
		/*{Description*/
		/*!Calls System1D<Type>(spuc,2,filename) to construct a system with
		 * spuc site per unit cell and 2 links per sites */
		/*}*/
		Chain(unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Chain()=0;

		Vector<double> compute_J(Vector<double> const& Jp);

	protected:
		void set_obs(int nobs);
		/*!Returns the neighbours of site i*/
		Matrix<int> get_neighbourg(unsigned int const& i) const;
		/*!Given N and m, save the best simulation in a text file for any n*/
		std::string extract_level_3();
		/*!Find the best range to compute the critcal exponents*/
		bool compute_critical_exponents(Vector<double> const& lrc, unsigned int& xi, unsigned int& xf, Vector<double>& p);

		virtual void energy_bound(){};
		void long_range_correlation_and_structure_factor();

	private:
		/*{Description*/
		/*!Do the fit of lrc over the range [xi,xf] then modify the parameters
		 * and check the quality of the fit*/
		/*}*/
		void do_fit(Vector<double> const& lrc, unsigned int const& xi, unsigned int const& xf, Vector<double>& p, double& R_squared, double& d_squared);
};

template<typename Type>
Chain<Type>::Chain(unsigned int const& spuc, std::string const& filename):
	System1D<Type>(spuc,2,filename)
{
	if(this->status_==2){
		if(!this->obs_.size()){
			this->set_nn_links(Vector<unsigned int>(1,1));
		}

		if(this->obs_[0].nlinks() != this->J_.size() && this->J_.size() != 0){
			Vector<double> tmp(this->J_);
			this->J_.set(this->obs_[0].nlinks());
			for(unsigned int i(0);i<this->J_.size();i++){ this->J_(i) = tmp(i%tmp.size()); }
		}

		if(this->J_.size()==this->obs_[0].nlinks()){
			std::string tmp("J");
			for(unsigned int i(0);i<this->spuc_;i++){ 
				tmp += ((this->J_(i)>0)?"+":"")+my::tostring(this->J_(i)); 
			}
			this->filename_.replace(this->filename_.find("Juniform"),8,tmp);
			this->path_.replace(this->path_.find("Juniform"),8,tmp);
		} else { std::cerr<<__PRETTY_FUNCTION__<<" : J_ has an incoherent size"<<std::endl; }
	}
}

template<typename Type>
Chain<Type>::~Chain() = default;

template<typename Type>
void Chain<Type>::set_obs(int nobs){
	if(nobs<0){ nobs = 2; }
	unsigned int nlinks;
	unsigned int nval;
	unsigned int m;
	if(nobs>0){/*bond energy*/
		nlinks = this->obs_[0].nlinks();
		nval = this->spuc_;
		this->obs_.push_back(Observable("Bond energy",1,nval,0));
		for(unsigned int i(0);i<nlinks;i++){
			this->obs_[0](i,2) = i%nval;
		}
	}
	if(nobs==2){/*long range correlation*/
		m = this->n_;
		nval = this->n_;
		nlinks = m*nval;
		this->obs_.push_back(Observable("Long range correlations",2,nval,nlinks));
		for(unsigned int i(0);i<m;i++){
			for(unsigned int j(0);j<nval;j++){
				this->obs_[2](i*nval+j,0) = i%this->n_;
				this->obs_[2](i*nval+j,1) = (i+j)%this->n_;
				this->obs_[2](i*nval+j,2) = j;
			}
		}
	}
}

template<typename Type>
Matrix<int> Chain<Type>::get_neighbourg(unsigned int const& i) const {
	Matrix<int> nb(this->z_,2);
	if(i!=this->n_-1){
		nb(0,0) = i+1; 
		nb(0,1) = 0;
	} else {
		nb(0,0) = 0;
		nb(0,1) = 1;
	}
	if(i!=0){
		nb(1,0) = i-1; 
		nb(1,1) = 0;
	} else {
		nb(1,0) = this->n_-1;
		nb(1,1) = 1;
	}
	return nb;
}

template<typename Type>
std::string Chain<Type>::extract_level_3(){
	double polymerization_strength;
	Vector<double> exponents;
	(*this->read_)>>polymerization_strength>>exponents;
	(*this->data_write_)<<this->N_<<" "<<this->m_<<" "<<this->bc_<<" "<<this->n_<<" "<<this->obs_[0][0]<<" "<<polymerization_strength<<" "<<exponents<<IOFiles::endl;

	return this->filename_;
}

template<typename Type>
bool Chain<Type>::compute_critical_exponents(Vector<double> const& lrc, unsigned int& xi, unsigned int& xf, Vector<double>& p){
	unsigned int dx(this->N_/this->m_);
	if(this->bc_ == 1 && this->n_>4*dx){
		double R_squared;
		double d_squared;
		if(this->N_ % 2){
			xi = (this->N_-this->m_)/2;//strange behaviour but needed for su9m3. maybe try something else
			xf = this->n_-xi-1;
			do_fit(lrc, xi, xf, p, R_squared, d_squared);

			/*!if the quality of the fit is bad, the next valid starting point is
			 * selected. this may not reduce the intrinsec quality of the fit, but
			 * it will reduce the error (rss)
			 * r_squared = 1-rss/tss*(x.size()-1)/(x.size()-exponents.size());
			 * rss /= (x.size()-exponents.size());
			 * it is actually needed only for su(5) m=1
			 */
			if(R_squared<0.9995){
				xi += dx;           //for 5 should start at (N/m-1)/2-1+N/m
				xf = this->n_-xi-1;//for 5 should stop at n-si-2
				do_fit(lrc, xi, xf, p, R_squared, d_squared);
			}
			return true;
		} else {
			xi = dx;
			xf = this->n_-dx-1;
			do{
				do_fit(lrc, xi, xf, p, R_squared, d_squared);

				if(R_squared > 0.999 && d_squared/this->m_ < 1.6e-8){ return true; }
				/*!normally all the fits have a R_square higher than 0.999. But
				 * if it is not the case, I just need a fit that is not too bad
				 * and I don't care about d*/
				if(R_squared < 0.999 && R_squared > 0.995 && d_squared < 1e-7){ return true; }
				xf -= dx;
				xi += dx;
			} while (xf-xi>2*dx);
			xi = dx;
			xf = lrc.size()-dx;
			do_fit(lrc, xi, xf, p, R_squared, d_squared);

			std::cerr<<__PRETTY_FUNCTION__<<" : No fit found for : SU("<<this->N_<<") m="<<this->m_<<" n="<<this->n_<<" (fitting range assumed)"<<std::endl;
			return false;
		}
	} else {
		xi = 1;
		xf = lrc.size()-2;
		p.set(4,0);
		return false;
	}
}

template<typename Type>
void Chain<Type>::do_fit(Vector<double> const& lrc, unsigned int const& xi, unsigned int const& xf, Vector<double>& p, double& R_squared, double& d_squared){
	double rss(0.0);
	double tss(0.0);
	double ym;
	p.set(4,2);
	p(1) -= 2.0/this->N_;
	Vector<double> x(xf-xi);
	Vector<double> y(xf-xi);
	for(unsigned int i(0);i<x.size();i++){
		x(i) = i+xi;
		y(i) = lrc(i+xi);
	}
	ym = y.mean();

	if(p.size()==3){
		auto func = [this](double x, const double* p){
			return p[0]*cos(2*M_PI*x*this->m_/this->N_)*(pow(x,-p[1])+pow(this->n_-x,-p[1]))+p[2]*(pow(x,-2.0)+pow(this->n_-x,-2.0));
		};
		Fit(x,y,p,func);
		for(unsigned int i(0);i<x.size();i++){
			rss += my::norm_squared(y(i)-func(x(i),p.ptr()));
			tss += my::norm_squared(ym-y(i));
		}
	}
	if(p.size()==4){
		auto func = [this](double x, const double* p){
			return p[0]*cos(2*M_PI*x*this->m_/this->N_)*(pow(x,-p[1])+pow(this->n_-x,-p[1]))+p[2]*(pow(x,-p[3])+pow(this->n_-x,-p[3]));
		};
		Fit(x,y,p,func);
		for(unsigned int i(0);i<x.size();i++){
			rss += my::norm_squared(y(i)-func(x(i),p.ptr()));
			tss += my::norm_squared(ym-y(i));
		}
	}

	R_squared = 1-rss/tss*(x.size()-1)/(x.size()-p.size());
	d_squared = rss/(x.size()-p.size());
}

template<typename Type>
void Chain<Type>::long_range_correlation_and_structure_factor(){
	if(this->obs_.size()>1){
		/*!long range correlations*/
		/*{*/
		IOFiles lr_corr_file(this->analyse_+this->path_+this->dir_+this->filename_+"-lr-c.dat",true);
		lr_corr_file<<"%j corr(i,j) dx conv(0|1) #conv mean(0|1)"<<IOFiles::endl;

		Vector<double> lr_corr(this->obs_[1].nval());
		for(unsigned int i(0);i<this->obs_[1].nval();i++){
			lr_corr_file<<i<<" "<<this->obs_[1][i]<<IOFiles::endl;
			lr_corr(i) = this->obs_[1][i].get_x();
		}

		unsigned int xi;
		unsigned int xf;
		Vector<double> exponents;
		bool fit(this->compute_critical_exponents(lr_corr,xi,xf,exponents));

		Gnuplot gplr(this->analyse_+this->path_+this->dir_,this->filename_+"-lr");
		gplr.range("x",this->N_/this->m_,this->n_-this->N_/this->m_);
		gplr.label("x","$\\|i-j\\|$","offset 0,0.5");
		gplr.label("y2","$<S_{\\alpha}^{\\alpha}(i)S_{\\alpha}^{\\alpha}(j)>-\\dfrac{m^2}{N}$","offset 1");
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
		gplr+="fit [" + my::tostring(xi) + ":" + my::tostring(xf) + "] f(x) '"+this->filename_+"-lr-c.dat' u 1:2 noerrors via p0,p1,p2,p3";
		gplr+="plot '"+this->filename_+"-lr-c.dat' u 1:2:3 w errorbars lt 1 lc 7 notitle,\\";
		gplr+="     f(x) lc 7 " + std::string(fit?"lw 0.5":"dt 2") + " t sprintf('$\\eta=%f$, $\\mu=%f$',p1,p3)";
		gplr.save_file();
		gplr.create_image(true,true);
		/*}*/
		/*!structure factor*/
		/*{*/
		unsigned int llr(this->obs_[1].nval());
		Vector<std::complex<double> > Ck(llr,0.0);
		std::complex<double> normalize(0.0);
		double dk(2.0*M_PI/llr);

		for(unsigned int k(0);k<llr;k++){
			for(unsigned int i(0);i<llr;i++){
				Ck(k) += std::polar(lr_corr(i),dk*k*i);
			}
			normalize += Ck(k);
		}
		Ck /= dk*normalize;

		IOFiles data_sf(this->analyse_+this->path_+this->dir_+this->filename_+"-sf.dat",true);
		for(unsigned int k(0);k<llr;k++){
			data_sf<<dk*k<<" "<<Ck(k).real()<<" "<<Ck(k).imag()<<IOFiles::endl;
		}

		Gnuplot gpsf(this->analyse_+this->path_+this->dir_,this->filename_+"-sf");
		gpsf+="set key bottom";
		gpsf.range("x","0","2*pi");
		switch(this->N_/this->m_){
			case 3: { gpsf+="set xtics ('0' 0,'$2\\pi/3$' 2.0*pi/3.0, '$4\\pi/3$' 4.0*pi/3.0,'$2\\pi$' 2.0*pi)"; } break;
			case 5: { gpsf+="set xtics ('0' 0,'$2\\pi/5$' 2.0*pi/5.0, '$4\\pi/5$' 4.0*pi/5.0, '$6\\pi/5$' 6.0*pi/5.0, '$8\\pi/5$' 8.0*pi/5.0, '$2\\pi$' 2.0*pi)"; } break;
			default:{ gpsf+="set xtics ('0' 0,'$\\pi/2$' pi/2.0,'$\\pi$' pi,'$3\\pi/2$' 3.0*pi/2.0,'$2\\pi$' 2.0*pi)"; } break;
		}
		gpsf.label("x","$k$","offset 0,0.5");
		gpsf.label("y2","$<S(k)>$");
		gpsf+="plot '"+this->filename_+"-sf.dat' u 1:2 lt 1 lc 6 t 'real',\\";
		gpsf+="     '"+this->filename_+"-sf.dat' u 1:3 lt 1 lc 7 t 'imag'";
		gpsf.save_file();
		gpsf.create_image(true,true);
		/*}*/

		//if(this->jd_write_){ this->jd_write_->write("critical exponents",exponents); }
	}
}
#endif
