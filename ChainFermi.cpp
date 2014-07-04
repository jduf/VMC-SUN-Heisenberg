#include"ChainFermi.hpp"

ChainFermi::ChainFermi(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc):
	System(ref,N,m,n,M,bc),
	Chain<double>(1,"chain-fermi")
{
	if(status_==1){
		init_fermionic();
		compute_T();

		system_info_.text("Spin ChainFermi, all the hopping parameters are real");
	}
}

/*{method needed for running*/
void ChainFermi::compute_T(){
	double t(1.0);
	T_.set(n_,n_,0);
	Matrix<int> nb;
	for(unsigned int i(0); i< n_; i++){
		nb = get_neighbourg(i);
		T_(i,nb(0,0)) = nb(0,1)*t;
	}
	T_ += T_.transpose();
}

void ChainFermi::create(){
	diagonalize_T();
	for(unsigned int c(0);c<N_;c++){
		if(!is_degenerate(c)){
			EVec_[c].set(n_,M_(c));
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = T_(i,j);
				}
			}
		}
	}
}
/*}*/

/*{method needed for checking*/
void ChainFermi::check(){
	BandStructure<double> bs(T_,Lx_,spuc_,bc_);
}
/*}*/
