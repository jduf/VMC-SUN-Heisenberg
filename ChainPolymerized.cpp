#include "ChainPolymerized.hpp"

ChainPolymerized::ChainPolymerized(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M, int const& bc, double delta):
	System(ref,N,m,n,M,bc),
	Chain<double>(n*m/N,"chain-polymerized"),
	delta_(delta)
{
	std::cout<<M_<<std::endl;
	init_fermionic();
	compute_T();

	filename_ += "-delta" + tostring(delta_);
	rst_.text("Spin chain, with different real hopping term.");
	rst_.text("For N colors and m particules per sites, every");
	rst_.text("N/m, there is a weaker bound, namely t-delta");
	rst_.text("instead of t+delta. (t=1,delta>0)");

	std::cout<<"chainpolymerized"<<std::endl;
}

ChainPolymerized::~ChainPolymerized(){}

void ChainPolymerized::create(unsigned int const& type){
	corr_.set(n_);
	if(type==2){ long_range_corr_.set(n_/3); }

	diagonalize_T('S');
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

void ChainPolymerized::compute_P(Matrix<double>& P){
	P.set(n_,n_,0);
	for(unsigned int i(0);i<N_/m_;i++){
		P(n_ -i-1,N_/m_-1-i) = bc_;
	}
	for(unsigned int i(0); i< n_-N_/m_; i++){
		P(i,i+N_/m_) = 1.0;
	}
}

void ChainPolymerized::check(){
	delta_=0.1;
	Matrix<double> P;
	compute_P(P);
	BandStructure<double> bs(T_,P);
	std::cout<<T_<<std::endl;
}

void ChainPolymerized::save(IOFiles& w) const{
	GenericSystem<double>::save(w);
	w("delta (t+-delta)",delta_);
}
