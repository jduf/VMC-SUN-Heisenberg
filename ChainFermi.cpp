#include"ChainFermi.hpp"

ChainFermi::ChainFermi(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc):
	System(ref,N,m,n,M,bc),
	Chain<double>(n,"chain-fermi")
{
	init_fermionic();
	compute_T();

	rst_.text("Spin ChainFermi, all the hopping parameters are real");
}

ChainFermi::~ChainFermi(){}

void ChainFermi::create(unsigned int const& type){
	std::cout<<"ChainFermi::create "<<" "<<type<<std::endl;
	T_.set(n_,n_,0);

	compute_T();
	diagonalize_T('S');
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
	compute_T();
	std::cout<<T_<<std::endl;
}
