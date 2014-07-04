#include "ChainPolymerized.hpp"

ChainPolymerized::ChainPolymerized(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc, double delta):
	System(ref,N,m,n,M,bc),
	Chain<double>(N/m,"chain-polymerized"),
	delta_(delta)
{
	if(status_==1){
		init_fermionic();
		compute_T();

		filename_ += "-delta" + tostring(delta_);
		system_info_.text("Spin chain, with different real hopping term.");
		system_info_.text("For N colors and m particules per sites, every");
		system_info_.text("N/m, there is a weaker bound, namely t-delta");
		system_info_.text("instead of t+delta. (t=1,delta>0)");
	}
}

/*{method needed for running*/
void ChainPolymerized::compute_T(){
	/*!If t<0, delta<0 otherwise no polymerization occurs
	 * If t>0, delta>0 otherwise no polymerization occurs */
	double t(1.0);
	T_.set(n_,n_,0);
	Matrix<int> nb;
	unsigned int a(n_/Lx_);
	for(unsigned int i(0); i < n_; i += a){
		for(unsigned int j(0); j<a; j++){
			nb = get_neighbourg(i+j);
			T_(i+j,nb(0,0)) = t+delta_;
		}
		nb = get_neighbourg(i+a-1);
		T_(i+a-1,nb(0,0)) = nb(0,1)*(t-delta_);
	}
	T_ += T_.transpose();
}

void ChainPolymerized::create(){
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);
	//if(type==2){ long_range_corr_.set(n_/3); }

	diagonalize_T();
	for(unsigned int c(0);c<N_;c++){
		if(!is_degenerate(c)){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = T_(i,j);
				}
			}
		}
	}
}

void ChainPolymerized::save(IOFiles& w) const{
	GenericSystem<double>::save(w);
	w("delta (t+-delta)",delta_);
}
/*}*/

/*{method needed for checking*/
void ChainPolymerized::check(){
	BandStructure<double> bs(T_,Lx_,spuc_,bc_);
}
/*}*/
