#ifndef DEF_CHAIN
#define DEF_CHAIN

#include "System1D.hpp"
#include "Fit.hpp"

template<typename Type>
class Chain: public System1D<Type>{
	public:
		/*{Description*/
		/*!Constructor of a system with spuc sites per unit cell. Calls the
		 * GenericSystem<Type>(2,filename), to construct a system with 2 links
		 * per sites */
		/*}*/
		Chain(unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Chain()=0;

	protected:
		/*!Returns the neighbours of site i*/
		Matrix<int> get_neighbourg(unsigned int i) const;
		/*!Find the best range to compute the critcal exponents*/
		void compute_critical_exponents(unsigned int& xi, unsigned int& xf, Vector<double>& exponents, Vector<double> const& lrc);
		std::string extract_level_3();
};

template<typename Type>
Chain<Type>::Chain(unsigned int const& spuc, std::string const& filename):
	System1D<Type>(spuc,2,filename)
{
	if(this->status_==1){ this->compute_links(); }
}

template<typename Type>
Chain<Type>::~Chain(){}

template<typename Type>
Matrix<int> Chain<Type>::get_neighbourg(unsigned int i) const {
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
void Chain<Type>::compute_critical_exponents(unsigned int& xi, unsigned int& xf, Vector<double>& exponents, Vector<double> const& lrc){
	double d(10);
	double bd(10);
	Vector<double> x;
	Vector<double> y;
	Vector<double> p(4,2);
	p(1) -= 2.0/this->N_;

	for(unsigned int s(1);s<lrc.size()/3;s++){
		x.set(lrc.size()-2*s+1);
		y.set(lrc.size()-2*s+1);
		for(unsigned int i(0);i<x.size();i++){
			x(i) = i+s;
			y(i) = lrc(i+s);
		}
		auto func = [this](double x, const double* p){ 
			return p[0]*cos(2*M_PI*x*this->m_/this->N_)*(pow(x,-p[1])+pow(this->n_-x,-p[1]))+p[2]*(pow(x,-p[3])+pow(this->n_-x,-p[3]));
		};
		Fit(x,y,p,func);
		d = norm_squared(p(1)-2.0+2.0/this->N_);
		if(d<bd){
			bd = d;
			xi = s;
			xf = this->n_-s; 
			exponents = p;
		}
	}
}

template<typename Type>
std::string Chain<Type>::extract_level_3(){
	double polymerization_strength;
	Vector<double> exponents;
	(*this->read_)>>this->E_>>polymerization_strength>>exponents;
	(*this->data_write_)<<this->N_<<" "<<this->m_<<" "<<this->bc_<<" "<<this->n_<<" "<<this->E_<<" "<<polymerization_strength<<" "<<exponents<<IOFiles::endl;

	return this->filename_;
}
#endif
