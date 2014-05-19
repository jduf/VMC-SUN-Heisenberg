#include"ChainFermi.hpp"

ChainFermi::ChainFermi(unsigned int N, unsigned int n, unsigned int m, int bc):
	Chain<double>(N,n,m,bc,"chain-fermi")
{
	rst_.text("Spin ChainFermi, all the hopping parameters are real");
}

ChainFermi::~ChainFermi(){}

bool ChainFermi::create(double x){
	compute_T();
	diagonalize_T('S');
	for(unsigned int spin(0);spin<N_;spin++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_;j++){
				EVec_(i+spin*n_,j) = T_(i,j);
			}
		}
	}
	return !degenerate_;
}

void ChainFermi::compute_T(){
	double t(-1.0);
	Matrix<int> nb;
	for(unsigned int i(0); i< n_; i++){
		nb = get_neighbourg(i);
		T_(i,nb(0,0)) = nb(0,1)*t;
	}
	T_ += T_.transpose();
}

void ChainFermi::compute_P(Matrix<double>& P){
	P.set(n_,n_);
	P(n_ -1,0) = bc_;
	for(unsigned int i(0); i< n_-1; i++){
		P(i,i+1) = 1.0;
	}
}

void ChainFermi::check(){
	bc_ = -1;
	compute_T();
	std::cout<<T_<<std::endl;
}
