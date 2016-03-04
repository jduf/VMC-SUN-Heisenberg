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
		/*!Plot the band structure E(p)*/
		void plot_band_structure();

	private:
		Matrix<Type> T_;	//!< translation operator along x-axis
		Vector<double> p_;	//!< eigenvalue of T
		Vector<double> e_;	//!< eigenvalue of H_

		/*!Returns the index of the site i in the unit cell*/
		unsigned int get_site_in_unit_cell(unsigned int const& i) const { return i%this->spuc_; }

		/*!Compute the translation operator*/
		void compute_T();
		/*!Diagonalize H_+T_ => compute the band structure E(p)*/
		bool full_diagonalization();
};
	
/*{constructors*/
template<typename Type>
System1D<Type>::System1D(unsigned int const& spuc, unsigned int const& z, std::string const& filename):
	GenericSystem<Type>(spuc,z,filename)
{
	if(this->N_%this->m_){ std::cerr<<"System1D : maybe problematric, m doesn't divide N, so check everywhere in the code where N/m appears"<<std::endl; }
	else { this->status_--; }
	//if(spuc%(this->N_/this->m_)){ std::cerr<<"System1D : problem in the definition of the unit cell"<<std::endl; }
	//else { this->status_--; } //don't understand why this line is important
}

template<typename Type>
System1D<Type>::~System1D() = default;
/*}*/

/*{protected methods*/
template<typename Type>
void System1D<Type>::plot_band_structure(){
	if(full_diagonalization()){
		IOFiles spectrum("spectrum.dat",true);
		for(unsigned int i(0);i<this->n_;i++){
			spectrum<<p_(i)<<" "<<e_(i)<<" "<<(i<this->M_(0))<<IOFiles::endl;
		}

		Gnuplot gp("./","spectrum");
		gp.range("x","-pi","pi");
		gp+="plot 'spectrum.dat' u ($3==1?$1:1/0):2 w p ps 1.5 lt 1 lc 4 t 'selected ev',\\";
		gp+="     'spectrum.dat' u ($3==0?$1:1/0):2 w p ps 1.5 lt 1 lc 7 notitle";
		gp.save_file();
		gp.create_image(true,false);
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : diagonalization failed, the band structure can't be computed"<<std::endl; }
}
/*}*/

/*{private methods*/
template<typename Type>
void System1D<Type>::compute_T(){
	T_.set(this->n_,this->n_,0);
	unsigned int tmp;
	unsigned int const L(this->n_/this->spuc_);
	for(unsigned int i(0); i<L-1; i++){
		tmp = this->spuc_*i;
		for(unsigned int k(0);k<this->spuc_;k++){
			T_(tmp+k,tmp+k+this->spuc_) = 1;
		}
	}
	tmp = this->n_-this->spuc_;
	for(unsigned int k(0);k<this->spuc_;k++){ T_(tmp+k,k) = this->bc_; }
}

template<typename Type>
bool System1D<Type>::full_diagonalization(){
	compute_T();
	Matrix<Type> M(this->H_);
	M += T_*Type(3.0);
	Vector<std::complex<double> > eval;
	Lapack<Type>(M,false,'G').eigensystem(eval,&this->evec_);

	for(unsigned int i(0);i<this->n_;i++){
		for(unsigned int j(i+1);j<this->n_;j++){
			if(my::are_equal(eval(i),eval(j),this->eq_prec_,this->eq_prec_)){
				std::cerr<<__PRETTY_FUNCTION__<<" : eigenvalue "<<i<<" and "<<j<<" degenerate"<<std::endl;
				return false;
			}
		}
	}
	Vector<unsigned int> index;
	e_.set(this->n_);
	for(unsigned int i(0);i<this->n_;i++){ e_(i) = this->projection(this->H_,i).real(); }
	e_.sort(std::less_equal<double>(),index);

	Matrix<std::complex<double> > evec_tmp(this->evec_);
	Vector<std::complex<double> > eval_tmp(eval);
	for(unsigned int i(0);i<this->n_;i++){
		for(unsigned int j(0);j<this->n_;j++){
			std::swap(this->evec_(i,j),evec_tmp(i,index(j)));
		}
		std::swap(eval(i),eval_tmp(index(i)));
	}

	p_.set(this->n_);
	for(unsigned int i(0);i<this->n_;i++){
		p_(i) = log(this->projection(T_,i)).imag();
	}
	return true;
}
/*}*/
#endif
