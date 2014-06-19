#include "ChainPolymerized.hpp"

ChainPolymerized::ChainPolymerized(unsigned int const& N, unsigned int const& n, unsigned int const& m, int const& bc, Vector<unsigned int> const& ref):
	System(N,n,m,bc,ref),
	Chain<double>("chain-polymerized")
{
	std::cout<<"chainpolymerized"<<N<<" "<<N_<<std::endl;
	rst_.text("Spin chain, with different real hopping term.");
	rst_.text("For N colors and m particules per sites, every");
	rst_.text("N/m, there is a weaker bound, namely t-delta");
	rst_.text("instead of t+delta. (t=1,delta>0)");
}

ChainPolymerized::~ChainPolymerized(){}

void ChainPolymerized::create(double const& delta, unsigned int const& type){
	delta_=delta;
	corr_.set(n_);
	if(type==2){ long_range_corr_.set(n_/3); }

	EVec_.set(n_*N_,M_,0);

	compute_T();
	diagonalize_T('S');
	for(unsigned int spin(0);spin<N_;spin++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_;j++){
				EVec_(i+spin*n_,j) = T_(i,j);
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
	for(unsigned int i(0); i < n_; i += a_){
		for(unsigned int j(0); j<a_-1; j++){
			nb = get_neighbourg(i+j);
			T_(i+j,nb(0,0)) = t+delta_;
		}
		nb = get_neighbourg(i+a_-1);
		T_(i+a_-1,nb(0,0)) = nb(0,1)*(t-delta_);
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

void ChainPolymerized::save(IOFiles& w) const{
	GenericSystem<double>::save(w);
	w("delta (t+-delta)",delta_);
}

void ChainPolymerized::check(){
	delta_=0.1;
	compute_T();
	Matrix<double> P;
	compute_P(P);
	BandStructure<double> bs(T_,P);
}
