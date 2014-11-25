#ifndef DEF_SYSTEM1D
#define DEF_SYSTEM1D

#include "GenericSystem.hpp"

/*{Description*/
/*!The main goal of this class is the selection of eigenvectors that will give
 * a minimal energy. The secondary one, is the visualization of the band
 * structure.
 *
 * For a ChainFermi and ChainPolymerized, the band structure is even and
 * each level is either non-degenerate or doubly degenerate. Therefore, for
 * a given energy level, there is no need to look for the eigenvector with
 * minimal impulsion as both of them will have the same absolute value.
 * !!!WRONG!!!
 *
 * When there is a degeneracy, the following happens : 
 *
 * + select all eigenvectors corresponding to an energy below E_F
 * + select the same eigenvector (as the same impulsion) for all colors
 *
 * It has been checked that this method works for 
 *
 * + SU(2) m=1
 * + SU(3) m=1
 * + SU(4) m=2 (works for n>=40)
 *
 * It has been observed that for those test, the energy found by this method is
 * the one found when the boundary condition are modified. The same could be
 * said about the structure factor altough it is less clear.
 */
/*}*/
template<typename Type>
class System1D: public GenericSystem<Type>{
	public:
		/*!Constructor*/
		System1D(unsigned int const& spuc, unsigned int const& z, std::string const& filename);
		/*!Destructor*/
		virtual ~System1D()=0;

	protected:
		Matrix<Type> H_;		//!< matrix used to get the band structure
		unsigned int const L_;	//!< number of unit cell along the x-axis
		Matrix<std::complex<double> > evec_;//!< eigenvector of H+T

		/*!Plot the band structure E(p)*/
		void plot_band_structure();
		/*!Create the selection of optimal eigenvectors*/
		void select_eigenvectors(unsigned int const& m);

	private:
		Matrix<Type> T_;	//!< translation operator along x-axis
		Vector<double> p_;	//!< eigenvalue of T
		Vector<double> e_;	//!< eigenvalue of H_

		/*!Compute the translation operator*/
		void compute_T();
		/*!Compute the band structure E(p)*/
		void compute_band_structure();
		/*!Evaluate the value of an operator O as <bra|O|ket>*/
		std::complex<double> projection(Matrix<Type> const& O, unsigned int const& idx);
};
	
/*{constructors*/
template<typename Type>
System1D<Type>::System1D(unsigned int const& spuc, unsigned int const& z, std::string const& filename):
	GenericSystem<Type>(spuc,z,filename),
	L_(this->n_/spuc)
{
	if(this->N_%this->m_){std::cout<<"System1D : maybe problematric, m doesn't devide N, so check everywhere in the code where N/m appears"<<std::endl;}
	this->status_--;
}

template<typename Type>
System1D<Type>::~System1D(){}
/*}*/

/*{protected methods*/
template<typename Type>
void System1D<Type>::select_eigenvectors(unsigned int const& m){
	if(!p_.size()){
		compute_T();
		compute_band_structure();
	}

	double n1(0);
	double n2(0);
	std::complex<double> tmp1;
	std::complex<double> tmp2;
	for(unsigned int i(0);i<this->n_;i++){
		tmp1 = evec_(i,m) + evec_(i,m-1);//k=k1+k2=0
		tmp2 = evec_(i,m) - evec_(i,m-1);//k=k1-k2=2k1
		evec_(i,m-1) = tmp1;
		evec_(i,m) = tmp2;
		n1 += norm_squared(tmp1);
		n2 += norm_squared(tmp2);
	}
	for(unsigned int i(0);i<this->n_;i++){
		evec_(i,m-1)/= sqrt(n1);
		evec_(i,m)  /= sqrt(n2);
	}
}

template<typename Type>
void System1D<Type>::plot_band_structure(){
	if(!p_.size()){
		compute_T();
		compute_band_structure();
	}

	IOFiles spectrum("spectrum.dat",true);
	for(unsigned int i(0);i<this->n_;i++){
		spectrum<<p_(i)<<" "<<e_(i)<<IOFiles::endl;
	}

	Gnuplot gp("./","spectrum");
	gp.xrange("-pi","pi");
	gp+="plot 'spectrum.dat' u 1:2 w p ps 1.5 lt 3 lc 7";
	gp.save_file();
}
/*}*/

/*{private methods*/
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
	Matrix<Type> M(H_);
	M += T_*Type(3.0);
	Vector<std::complex<double> > eval;
	Lapack<Type>(M,true,'G').eigensystem(eval,&evec_);

	Vector<unsigned int> index;
	e_.set(this->n_);
	for(unsigned int i(0);i<this->n_;i++){ e_(i) = projection(H_,i).real(); }
	e_.sort(std::less_equal<double>(),index);

	Matrix<std::complex<double> > evec_tmp(evec_);
	Vector<std::complex<double> > eval_tmp(eval);
	for(unsigned int i(0);i<this->n_;i++){
		for(unsigned int j(0);j<this->n_;j++){
			std::swap(evec_(i,j),evec_tmp(i,index(j)));
		}
		std::swap(eval(i),eval_tmp(index(i)));
	}

	p_.set(this->n_);
	for(unsigned int i(0);i<this->n_;i++){
		p_(i) = log(projection(T_,i)).imag();
	}

	for(unsigned int i(0);i<this->n_;i++){
		for(unsigned int j(i+1);j<this->n_;j++){
			if(are_equal(eval(i),eval(j),1e-10,1e-10)){
				this->degenerate_ = true;
				std::cout<<"H+T eigenvalue degenerate"<<std::endl;
				i=j=this->n_;
			}
		}
	}
	if(this->degenerate_){ std::cerr<<"void System1D<Type>::compute_band_structure() : degenerate"<<std::endl; }

	std::cout<<"val"<<std::endl;
	std::cout<<eval.chop()<<std::endl;
	std::cout<<"p"<<std::endl;
	std::cout<<p_<<std::endl;
	std::cout<<"e"<<std::endl;
	std::cout<<e_<<std::endl;
	////std::cout<<T_<<std::endl;
	//std::cout<<H_<<std::endl;
	//std::cout<<"vec"<<std::endl;
	//std::cout<<evec_.chop()<<std::endl;
	

}

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
