#include "KagomeFermi.hpp"

template<>
void KagomeFermi<double>::create(){
	std::cout<<"double"<<std::endl;
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);

	compute_T();
	diagonalize_T();

	if(!degenerate_){
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
void KagomeFermi<std::complex<double> >::create(){
	std::cout<<"complex"<<std::endl;
	E_.set(50,5,false);
	corr_.set(links_.row(),50,5,false);

	compute_T();
	BandStructure<std::complex<double> > bs(T_,Lx_,Ly_,spuc_,bc_);
	bs.compute_band_structure();
	bs.save();
	Matrix<std::complex<double> > evec(bs.get_evec());
	Vector<double> E(bs.get_E());
	Vector<double> px(bs.get_px());
	Vector<double> py(bs.get_py());
	Vector<unsigned int> index;
	E.sort(std::less_equal<double>(),index);

	//bs.diagonalize_everything(T_,eval);
	//std::cout<<eval.chop()<<std::endl;
	//unsigned int a(0);
	//unsigned int b(0);
	//while(a != n_){
	//do{b++;}
	//while(b != n_ && are_equal(eval_(b),eval_(b-1),1e-14));
	//bs.diagonalize_subspace_Tx(a,b,T_,eval_);
	//a=b;
	//}
	//bs.save();
	
	for(unsigned int i(0);i<n_;i++){
		std::cout<<E(i)<<(i==M_(0)?"|":" ");
	}
	std::cout<<std::endl;

	for(unsigned int c(0);c<N_;c++){
		EVec_[c].set(n_,M_(c));
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_(c)-1;j++){
				EVec_[c](i,j) = evec(i,index(j));
			}
			EVec_[c](i,M_(c)-1) = evec(i,index(M_(c)-1+c));
		}
		double Px(px(index(M_(c)-1)));
		double Py(py(index(M_(c)-1)));
		double e(  E(index(M_(c)-1)));
		for(unsigned int j(0);j<M_(c)-1;j++){
			Px += (are_equal(std::abs(px(index(j))),M_PI,1e-8)?0:px(index(j)));
			Py += (are_equal(std::abs(py(index(j))),M_PI,1e-8)?0:px(index(j)));
			e += E(index(j));
		}
		std::cout<<e<<" "<<Px<<" "<<Py<<std::endl;
	}
	T_.set();
}
