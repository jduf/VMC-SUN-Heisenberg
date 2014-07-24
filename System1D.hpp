#ifndef DEF_SYSTEM1D
#define DEF_SYSTEM1D

#include "GenericSystem.hpp"

template<typename Type>
class System1D: public GenericSystem<Type>{
	public:
		/*!Constructor*/
		System1D(unsigned int const& spuc, unsigned int const& z, std::string const& filename);
		/*!Destructor*/
		virtual ~System1D()=0;

	protected:
		Matrix<Type> H_;			//!< matrix used to get the band structure
		unsigned int const L_;		//!< number of unit cell along the x-axis
		Vector<unsigned int>* select_;
		Matrix<std::complex<double> > evec_;//!< eigenvector of H+T

		void compute_T();
		void compute_band_structure();
		void plot_band_structure();
		void select_eigenvectors();

	private:
		Matrix<Type> T_;	//!< translation operator along x-axis
		Vector<double> e_;	//!< eigenvalue of the Hamiltonian T_
		Vector<double> p_;	//!< eigenvalue of the translation along x
		unsigned int sel_[2];

		/*!Evaluate the value of an operator O as <bra|O|ket>*/
		std::complex<double> projection(Matrix<Type> const& O, unsigned int const& idx);
};
	
/*{constructors*/
template<typename Type>
System1D<Type>::System1D(unsigned int const& spuc, unsigned int const& z, std::string const& filename):
	GenericSystem<Type>(spuc,z,filename),
	L_(this->n_/spuc),
	select_(new Vector<unsigned int>[this->N_])
{
	sel_[0]= 0;
	sel_[1]= 0;
}

template<typename Type>
System1D<Type>::~System1D(){
	if(select_){ delete[] select_;}
}
/*}*/

/*{protected methods*/
template<typename Type>
void System1D<Type>::compute_T(){
	T_.set(this->n_,this->n_,0);
	unsigned int tmp;
	double t(1);
	for(unsigned int i(0); i<L_-1; i++){ 
		tmp = this->spuc_*i;
		for(unsigned int k(0);k<this->spuc_;k++){
			T_(tmp+k,tmp+k+this->spuc_) = t; 
		}
	}
	tmp = this->spuc_*(L_-1);
	for(unsigned int k(0);k<this->spuc_;k++){ T_(tmp+k,k) = this->bc_*t; }
}

template<typename Type>
void System1D<Type>::compute_band_structure(){
	e_.set(this->n_);
	Rand rnd(1e4);
	Matrix<Type> M;
	bool degenerate;
	unsigned int iter(0);
	do{
		/*may be optimized : */
		M  = H_;
		M += T_*Type(rnd.get()*rnd.get(1e2));

		Vector<std::complex<double> > eval;
		Lapack<Type>(M,true,'G').eigensystem(eval,&evec_);
		for(unsigned int i(0);i<this->n_;i++){
			e_(i) = projection(H_,i).real();
		}

		degenerate = false;
		for(unsigned int i(0);i<this->n_;i++){
			for(unsigned int j(i+1);j<this->n_;j++){
				if(are_equal(eval(i),eval(j),1e-10)){
					degenerate = true;
				}
			}
		}
	} while ( degenerate && iter++<1e2);
	if( iter>1e2 ){ std::cerr<<"void System1D<Type>::compute_band_structure() : degenerate"<<std::endl; }
	else {
		Vector<unsigned int> index;
		Matrix<std::complex<double> > evec_tmp(evec_);
		e_.sort(std::less_equal<double>(),index);
		p_.set(this->n_);
		for(unsigned int i(0);i<this->n_;i++){
			for(unsigned int j(0);j<this->n_;j++){
				std::swap(evec_(i,j),evec_tmp(i,index(j)));
			}
		}
		for(unsigned int i(0);i<this->n_;i++){
			p_(i) = log(projection(T_,i)).imag();
		}
	}
}

template<typename Type>
void System1D<Type>::plot_band_structure(){
	IOFiles spectrum("spectrum.dat",true);
	for(unsigned int i(0);i<this->n_;i++){
		spectrum<<p_(i)<<" "<<e_(i)<<IOFiles::endl;
	}

	Gnuplot gp("./","spectrum");
	gp.xrange("-pi","pi");
	gp+="plot 'spectrum.dat' u 1:2";
	gp.save_file();
}

template<typename Type>
void System1D<Type>::select_eigenvectors(){
	//std::cout<<e_<<std::endl;
	unsigned int iter(0);
	double P(0.0);
	for(unsigned int c(0);c<this->N_;c++){
		unsigned int a(this->M_(c)-1);
		unsigned int b(this->M_(c)-1);
		//std::cout<<a<<" "<<b<<std::endl;
		do{b++;} while (b+1<this->n_ && are_equal(e_(b),e_(b-1)));
		//std::cout<<a<<" "<<b<<std::endl;
		if(b!=this->M_(c)){ while(a>0 && are_equal(e_(a-1),e_(a))){a--;} }
		//std::cout<<a<<" "<<b<<std::endl;
		//std::cout<<e_.range(a,b)<<std::endl;
		Vector<unsigned int> cnk;
		Combination cbn;
		cbn.set(this->M_(c)-a,b-a,cnk);
		select_[c].set(this->M_(c));
		for(unsigned int i(0);i<a;i++){ select_[c](i) = i; }
		//do{ 
		//for(unsigned int i(0);i+a<M_(c);i++){ select_[c](a+i) = a+cnk(i); }
		//double P(0.0);
		//double  e(0.0);
		//for(unsigned int i(0);i<M_(c);i++){
		////Px+= (are_equal(std::abs(p_(a+comb(j))),M_PI,1e-8)?0:p_(a+comb(j)));
		//e += e_(select_[c](i));
		//}
		//} while (cbn.next());
		while(iter++<sel_[c] && cbn.next());
		for(unsigned int i(0);i+a<this->M_(c);i++){ select_[c](a+i) = a+cnk(i); }
		double  e(0.0);
		for(unsigned int i(0);i<this->M_(c);i++){
			P+= are_equal(std::abs(p_(select_[c](i))),M_PI,1e-12)?0:p_(select_[c](i));
			e +=  e_(select_[c](i));
		}
	}
}
/*}*/

/*{private methods*/
template<typename Type>
std::complex<double> System1D<Type>::projection(Matrix<Type> const& O, unsigned int const& idx){
	std::complex<double> tmp;
	std::complex<double> out(0.0);
	for(unsigned int i(0);i<O.row();i++){
		tmp = 0.0;
		for(unsigned int j(0);j<O.col();j++){
			tmp += O(i,j)*evec_(j,idx);
		}
		out += std::conj(evec_(i,idx))*tmp;
	}
	return out;
}
/*}*/
#endif
