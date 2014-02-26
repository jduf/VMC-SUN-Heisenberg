#include"ChainPolymerized.hpp"

ChainPolymerized::ChainPolymerized(unsigned int N, unsigned int n, unsigned int m):
	Chain<double>(N,n,m,"chain-polymerized")
{
	rst_.text("Spin chain, with different hopping term on the odd and even sites");
}

ChainPolymerized::~ChainPolymerized(){}

void ChainPolymerized::create(double delta){
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
}

void ChainPolymerized::compute_T(){
	Vector<unsigned int> neighbourg;
	double t(1.0);
	for(unsigned int i(0); i < n_; i += a_){
		for(unsigned int j(0); j<a_-1; j++){
			neighbourg = get_neighbourg(i+j);
			T_(i+j,neighbourg(0)) = t+delta_;
		}
		neighbourg = get_neighbourg(i+a_-1);
		T_(i+a_-1,neighbourg(0)) = t-delta_;
	}
	T_(n_-1,0) = bc_*(t-delta_);
	T_ += T_.transpose();
}

void ChainPolymerized::compute_P(Matrix<double>& P){
	P.set(n_,n_);
	P(n_ -1,0) = bc_;
	for(unsigned int i(0); i< n_-1; i++){
		P(i,i+1) = 1.0;
	}
}

void ChainPolymerized::save(Write& w){
	GenericSystem<double>::save(w);
	w("delta (t+-delta)",delta_);
}
