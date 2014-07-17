#include"ChainFermi.hpp"

template<>
void ChainFermi<double>::create(){
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);
	//if(type==2){ long_range_corr_.set(n_/3); }

	compute_T();
	diagonalize_T();

	if(degenerate_){
		for(unsigned int c(0);c<N_;c++){
			EVec_[c].set(n_,M_(c));
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = T_(i,j);
				}
			}
		}
	}
}

template<>
void ChainFermi<std::complex<double> >::create(){
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);
	//if(type==2){ long_range_corr_.set(n_/3); }

	compute_T();
	BandStructure<std::complex<double> > bs(T_,Lx_,spuc_,bc_);
	diagonalize_T();

	unsigned int a(0);
	unsigned int b(0);
	while(a != n_){
		b = a;
		do{b++;}
		while(b != n_ && are_equal(eval_(b),eval_(b-1),1e-14));
		bs.dostuff(a,b,T_,eval_);
		a=b;
	}
	bs.save();

	Vector<double> P(bs.check(T_));
	for(unsigned int c(0);c<N_;c++){
		double p(0);
		double e(0);
		if(std::abs(eval_(M_(c)) - eval_(M_(c)-1))<1e-12){
			std::cerr<<"Fermi level degenerate"<<std::endl;
			//degenerate_=true;
		}
		EVec_[c].set(n_,M_(c));
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_(c)-1;j++){
				EVec_[c](i,j) = T_(i,j);
				p+=P(i);
				e+=eval_(i);
			}
			EVec_[c](i,M_(c)-1) = T_(i,M_(c)-1+c);
			p+=                     P(M_(c)-1+c);
			e+=                 eval_(M_(c)-1+c);
		}
		std::cout<<c<<" "<<e/double(n_)<<" "<<p/n_<<std::endl;
	}

	//std::cout<<P<<std::endl;
	//std::cout<<eval_<<std::endl;
}
