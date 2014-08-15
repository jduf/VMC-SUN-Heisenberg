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
 *
 * When there is a degeneracy, the following happens : 
 *
 * + select all eigenvectors corresponding to an energy below E_F
 * + select the same eigenvector (has the same impulsion) for all colors
 *
 * It has been checked that this method works for 
 *
 * + SU(2) m=1
 * + SU(3) m=1
 * + SU(4) m=2 (works for n>=40)
 *
 * It has been observed that for those test, the energy found by this method is
 * the one found when the boundary condition are modified. The same could be
 * said about the structure factor altough it is less clear
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
		Matrix<Type> H_;			//!< matrix used to get the band structure
		unsigned int const L_;		//!< number of unit cell along the x-axis
		Vector<unsigned int>* select_;//!< obtimal selection of eigenvectors
		Matrix<std::complex<double> > evec_;//!< eigenvector of H+T

		/*!Plot the band structure E(p)*/
		void plot_band_structure();
		/*!Create the selection of optimal eigenvectors*/
		void select_eigenvectors();

	private:
		Matrix<Type> T_;	//!< translation operator along x-axis
		Vector<double> e_;	//!< eigenvalue of the Hamiltonian T_
		Vector<double> p_;	//!< eigenvalue of the translation along x

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
	L_(this->n_/spuc),
	select_(new Vector<unsigned int>[this->N_])
{
	if(this->N_%this->m_){std::cout<<"System1D : maybe problematric, m doesn't devide N, so check everywhere in the code where N/m appears"<<std::endl;}
	this->status_--;
}

template<typename Type>
System1D<Type>::~System1D(){
	if(select_){ delete[] select_;}
}
/*}*/

/*{protected methods*/
template<typename Type>
void System1D<Type>::select_eigenvectors(){
	if(!p_.size()){
		compute_T();
		compute_band_structure();
	}

	for(unsigned int c(0);c<this->N_;c++){
		select_[c].set(this->M_(c));
		for(unsigned int i(0);i<this->M_(c);i++){ select_[c](i) = i; }
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
	gp+="plot 'spectrum.dat' u 1:2";
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

	Vector<unsigned int> index;
	e_.set(this->n_);
	for(unsigned int i(0);i<this->n_;i++){ e_(i) = projection(H_,i).real(); }
	e_.sort(std::less_equal<double>(),index);

	Matrix<std::complex<double> > evec_tmp(evec_);
	for(unsigned int i(0);i<this->n_;i++){
		for(unsigned int j(0);j<this->n_;j++){
			std::swap(evec_(i,j),evec_tmp(i,index(j)));
		}
	}
	
	p_.set(this->n_);
	for(unsigned int i(0);i<this->n_;i++){
		p_(i) = log(projection(T_,i)).imag();
	}
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
