#ifndef DEF_SYSTEM1D
#define DEF_SYSTEM1D

#include "GenericSystem.hpp"

/*{Description*/
/*!This class allow the diagonalization of the trial hamiltonian and the
 * visualization of the band structure.
 *
 * The band structure is computed as follows :
 *
 * + diagonalize H+3T
 * + use the eigenvectors to compute e,kx
 * + a linear combination of the degenerate eigenvectors |E_F,+>,|E_F,->
 * such that the new ones are |0>=|E_F,+>+|E_F,-> and |k>=|E_F,+>-|E_F,->
 */
/*}*/
template<typename Type>
class System1D: public GenericSystem<Type>{
	public:
		/*!Constructor*/
		System1D(unsigned int const& spuc, unsigned int const& z, std::string const& filename);
		/*!Default destructor*/
		virtual ~System1D()=0;

	protected:
		Matrix<Type> H_;		//!< matrix used to get the band structure
		unsigned int const L_;	//!< number of unit cell along the x-axis
		Matrix<std::complex<double> > evec_;//!< eigenvector of H+T

		/*!Plot the band structure E(p)*/
		void plot_band_structure();
		/*!Create the selection of optimal eigenvectors*/
		void select_eigenvectors();

		void diagonalize(bool simple);

	private:
		Matrix<Type> T_;	//!< translation operator along x-axis
		Vector<double> p_;	//!< eigenvalue of T
		Vector<double> e_;	//!< eigenvalue of H_

		/*!Compute the translation operator*/
		void compute_T();
		/*!Diagonalize H_*/
		bool simple_diagonalization();
		/*!Diagonalize H_+T_ => compute the band structure E(p)*/
		bool full_diagonalization();
		/*!Evaluate the value of an operator O as <bra|O|ket>*/
		std::complex<double> projection(Matrix<Type> const& O, unsigned int const& idx);
};
	
/*{constructors*/
template<typename Type>
System1D<Type>::System1D(unsigned int const& spuc, unsigned int const& z, std::string const& filename):
	GenericSystem<Type>(spuc,z,filename),
	L_(this->n_/spuc)
{
	if(this->N_%this->m_){ std::cout<<"System1D : maybe problematric, m doesn't devide N, so check everywhere in the code where N/m appears"<<std::endl; }
	this->status_--;
}

template<typename Type>
System1D<Type>::~System1D() = default;
/*}*/

/*{protected methods*/
template<typename Type>
void System1D<Type>::diagonalize(bool simple){
	if(simple){ if(simple_diagonalization()){ this->status_--; } }
	else { if(full_diagonalization()){ this->status_--; } }
}

template<typename Type>
void System1D<Type>::plot_band_structure(){
	if(full_diagonalization()){
		IOFiles spectrum("spectrum.dat",true);
		for(unsigned int i(0);i<this->n_;i++){
			spectrum<<p_(i)<<" "<<e_(i)<<IOFiles::endl;
		}

		Gnuplot gp("./","spectrum");
		gp.range("x","-pi","pi");
		gp+="plot 'spectrum.dat' u 1:2 w p ps 1.5 lt 3 lc 7";
		gp.save_file();
	} else {
		std::cerr<<__PRETTY_FUNCTION__<<" : diagonalization failed, the band structure can't be computed"<<std::endl;
	}
}
/*}*/

/*{private methods*/
template<typename Type>
void System1D<Type>::compute_T(){
	T_.set(this->n_,this->n_,0);
	unsigned int tmp;
	for(unsigned int i(0); i<L_-1; i++){ 
		tmp = this->spuc_*i;
		for(unsigned int k(0);k<this->spuc_;k++){
			T_(tmp+k,tmp+k+this->spuc_) = 1; 
		}
	}
	tmp = this->n_-this->spuc_;
	for(unsigned int k(0);k<this->spuc_;k++){ T_(tmp+k,k) = this->bc_; }
}

template<typename Type>
bool System1D<Type>::simple_diagonalization(){
	Vector<double> eval;
	Lapack<Type>(H_,false,(this->ref_(1)==1?'S':'H')).eigensystem(eval,true);
	for(unsigned int c(0);c<this->N_;c++){
		if(my::are_equal(eval(this->M_(c)),eval(this->M_(c)-1),1e-12)){
			std::cerr<<__PRETTY_FUNCTION__<<" : degenerate at the Fermi level"<<std::endl;
			return false;
		}
	}
	return true;
}

template<typename Type>
bool System1D<Type>::full_diagonalization(){
	compute_T();
	Matrix<Type> M(H_);
	M += T_*Type(3.0);
	Vector<std::complex<double> > eval;
	Lapack<Type>(M,true,'G').eigensystem(eval,&evec_);

	for(unsigned int i(0);i<this->n_;i++){
		for(unsigned int j(i+1);j<this->n_;j++){
			if(my::are_equal(eval(i),eval(j),1e-10,1e-10)){
				std::cerr<<__PRETTY_FUNCTION__<<" : eigenvalue "<<i<<" and "<<j<<" degenerate"<<std::endl;
				return false;
			}
		}
	}
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
	return true;
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
