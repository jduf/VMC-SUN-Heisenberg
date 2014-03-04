#include"ChainPolymerized.hpp"

ChainPolymerized::ChainPolymerized(unsigned int N, unsigned int n, unsigned int m, int bc):
	Chain<double>(N,n,m,bc,"chain-polymerized")
{
	rst_.text("Spin chain, with different hopping term on the odd and even sites");
}

ChainPolymerized::~ChainPolymerized(){}

unsigned int ChainPolymerized::create(double delta){
	filename_ += "-delta" + tostring(delta);
	delta_=delta;

	compute_T();
	diagonalize_T('S');
	for(unsigned int spin(0);spin<N_;spin++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_;j++){
				EVec_(i+spin*n_,j) = T_(i,j);
			}
		}
	}
	if(degenerate_){ return 0; }
	else { return 1; }/*1st step successful*/
}

void ChainPolymerized::compute_T(){
	/*!If t<0, delta<0 otherwise no polymerization occurs
	 * If t>0, delta>0 otherwise no polymerization occurs */
	double t(1.0);
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
	P.set(n_,n_);
	P(n_ -1,0) = bc_;
	for(unsigned int i(0); i< n_-1; i++){
		P(i,i+1) = 1.0;
	}
}

void ChainPolymerized::save(Write& w) const{
	GenericSystem<double>::save(w);
	w("delta (t+-delta)",delta_);
}

void ChainPolymerized::check(){
	delta_=0.1;
	compute_T();
	std::cout<<T_<<std::endl;
	std::cout<<sts_<<std::endl;
}
