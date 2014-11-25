#ifndef DEF_SYSTEM1D
#define DEF_SYSTEM1D

#include "GenericSystem.hpp"

/*{Description*/
/*!The main goal of this class is the selection of eigenvectors that will give
 * a minimal energy in agreement with periodic boundary condition. The
 * secondary one, is the visualization of the band structure.
 *
 * 1D chain with all hopping term having the same amplitude. The band structure
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
 * eigenvector unequivocal and there is no need to worry further. But when
 * delta==0, ChainFermi and ChainPolymerized are equivalent and suffers from
 * the same problem, the selection of eigenvector is equivocal because at the
 * "Fermi" level, the energies are degenerate. The only case when this can be
 * avoided is when (spuc && n) are even and n/spuc is odd because the band
 * structure (1) selects unequivocally the good eigenvectors. For all other
 * cases, the following method should be applied : 
 *
 * + diagonalize H+3T
 * + use the eigenvectors to compute e,kx
 * + make a linear combination of the degenerate eigenvectors |E_F,+>,|E_F,->
 * such that the new ones are |0>=|E_F,+>+|E_F,-> and |k>=|E_F,+>-|E_F,->
 * + |0> should be real an |k> complex
 * + select |0> to complete the selection of eigenvectors
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

	for(unsigned int i(0);i<this->n_;i++){
		for(unsigned int j(i+1);j<this->n_;j++){
			if(are_equal(eval(i),eval(j),1e-10,1e-10)){
				this->degenerate_ = true;
				std::cout<<"H+T eigenvalue degenerate"<<std::endl;
				i=j=this->n_;
			}
		}
	}
	if(this->degenerate_){
		std::cerr<<"void System1D<Type>::compute_band_structure() : degenerate"<<std::endl; 
	} else {
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
