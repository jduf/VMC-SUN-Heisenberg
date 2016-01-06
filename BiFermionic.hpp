#ifndef DEF_BIFERMIONIC
#define DEF_BIFERMIONIC

#include "System.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class BiFermionic : public virtual System{
	public:
		/*!Constructor*/
		BiFermionic(Fermionic<Type> const& F0, Fermionic<Type> const& F1);
		/*!Copy constructor*/
		BiFermionic(BiFermionic<Type> const& BF);
		/*!Constructor that reads from file*/
		BiFermionic(IOFiles& r);
		/*!Destructor*/
		virtual ~BiFermionic();
		/*{Forbidden*/
		BiFermionic(BiFermionic<Type>&&) = delete;
		BiFermionic<Type>& operator=(BiFermionic<Type>) = delete;
		/*}*/

		BiFermionic<Type> const& get_bifermionic() const { return (*this); }

	protected:
		/*!Constructor*/
		BiFermionic();

		bool same_wf_[2];		//!< true if the same wavefunction is used for all colors
		Matrix<Type>* EVec_[2]; //!< eigenvectors matrix (transfer matrix)
};

/*constructors and destructor and initialization*/
/*{*/
template<typename Type>
BiFermionic<Type>::BiFermionic(Fermionic<Type> const& F0, Fermionic<Type> const& F1):
	System(F1),
	same_wf_{F0.use_same_wf(),F1.use_same_wf()},
	EVec_{new Matrix<Type>[N_],new Matrix<Type>[N_]}
{
	if(F0.get_EVec() || F1.get_EVec()){
		if(same_wf_[0]){ for(unsigned int c(0);c<N_;c++){ EVec_[0][c] = F0.get_EVec(0); } }
		else { for(unsigned int c(0);c<N_;c++){ EVec_[0][c] = F0.get_EVec(c); } }
		if(same_wf_[1]){ for(unsigned int c(0);c<N_;c++){ EVec_[1][c] = F1.get_EVec(0); } }
		else { for(unsigned int c(0);c<N_;c++){ EVec_[1][c] = F1.get_EVec(c); } }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : EVec_ == NULL"<<std::endl; }
}

template<typename Type>
BiFermionic<Type>::BiFermionic(BiFermionic<Type> const& BF):
	System(BF),
	same_wf_{BF.same_wf_[0],BF.same_wf_[1]},
	EVec_{new Matrix<Type>[N_],new Matrix<Type>[N_]}
{
	if(BF.EVec_){
		for(unsigned int i(0);i<2;i++){
			if(same_wf_[i]){ for(unsigned int c(0);c<N_;c++){ EVec_[i][c] = BF.EVec_[i][0]; } }
			else { for(unsigned int c(0);c<N_;c++){ EVec_[i][c] = BF.EVec_[i][c]; } }
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : EVec_ == NULL"<<std::endl; }
}

template<typename Type>
BiFermionic<Type>::BiFermionic(IOFiles& r):
	System(r),
	same_wf_{r.read<bool>(),r.read<bool>()},
	EVec_{new Matrix<Type>[N_],new Matrix<Type>[N_]}
{
	if(N_){
		for(unsigned int i(0);i<2;i++){
			r>>EVec_[i][0];
			if(same_wf_[i]){ for(unsigned int c(1);c<N_;c++){ EVec_[i][c] = EVec_[i][0]; } }
			else { for(unsigned int c(1);c<N_;c++){ r>>EVec_[i][c]; } }
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : N_ == 0"<<std::endl; }
}

template<typename Type>
BiFermionic<Type>::BiFermionic():
	same_wf_{true,true},
	EVec_{NULL,NULL}
{}

template<typename Type>
BiFermionic<Type>::~BiFermionic(){
	if(EVec_[0]){ delete[] EVec_[0]; }
	if(EVec_[1]){ delete[] EVec_[1]; }
}
/*}*/
#endif
