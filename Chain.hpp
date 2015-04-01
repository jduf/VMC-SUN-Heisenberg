#ifndef DEF_CHAIN
#define DEF_CHAIN

#include "System1D.hpp"
#include "Fit.hpp"

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

	protected:
		/*!Returns the neighbours of site i*/
		Matrix<int> get_neighbourg(unsigned int const& i) const;
		/*!Given N and m, save the best simulation in a text file for any n*/
		std::string extract_level_3();
		/*!Find the best range to compute the critcal exponents*/
		bool compute_critical_exponents(Vector<double> const& lrc, unsigned int& xi, unsigned int& xf, Vector<double>& p);

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
		this->compute_links(); 
	}
}

template<typename Type>
Chain<Type>::~Chain(){}

template<typename Type>
Matrix<int> Chain<Type>::get_neighbourg(unsigned int const& i) const {
	Matrix<int> nb(this->z_,2,1);
	if( i != this->n_-1){ nb(0,0) = i+1;}
	else {
		nb(0,0) = 0;
		nb(0,1) = this->bc_;
	}
	if( i != 0){ nb(1,0) = i-1;}
	else {
		nb(1,0) = this->n_-1;
		nb(1,1) = this->bc_;
	}
	return nb;
}

template<typename Type>
std::string Chain<Type>::extract_level_3(){
	double polymerization_strength;
	Vector<double> exponents;
	(*this->read_)>>this->E_>>polymerization_strength>>exponents;
	(*this->data_write_)<<this->N_<<" "<<this->m_<<" "<<this->bc_<<" "<<this->n_<<" "<<this->E_<<" "<<polymerization_strength<<" "<<exponents<<IOFiles::endl;

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

			std::cerr<<"void Chain<Type>::compute_critical_exponents(unsigned int& xi,"<<std::endl;
			std::cerr<<"unsigned int& xf, Vector<double>& exponents, Vector<double> lrc) :"<<std::endl;
			std::cerr<< "    No fit found for : SU("<<this->N_<<") m="<<this->m_<<" n="<<this->n_;
			std::cerr<< " (fitting range assumed)"<<std::endl;
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
			rss += norm_squared(y(i)-func(x(i),p.ptr()));
			tss += norm_squared(ym-y(i));
		}
	}
	if(p.size()==4){
		auto func = [this](double x, const double* p){ 
			return p[0]*cos(2*M_PI*x*this->m_/this->N_)*(pow(x,-p[1])+pow(this->n_-x,-p[1]))+p[2]*(pow(x,-p[3])+pow(this->n_-x,-p[3]));
		};
		Fit(x,y,p,func); 
		for(unsigned int i(0);i<x.size();i++){
			rss += norm_squared(y(i)-func(x(i),p.ptr()));
			tss += norm_squared(ym-y(i));
		}
	}

	R_squared = 1-rss/tss*(x.size()-1)/(x.size()-p.size());
	d_squared = rss/(x.size()-p.size());
}
#endif
